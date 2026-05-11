import os
import glob
import cv2
import numpy as np
import pandas as pd
import torch
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
from collections import defaultdict
import re
import random
import xml.etree.ElementTree as ET
from PIL import Image
import segmentation_models_pytorch as smp
import albumentations as A
from albumentations.pytorch import ToTensorV2

IMAGE_CONFIG = {
    'img_height': 1024,
    'img_width': 1024,
    'input_channels_config': 3,
    'pixel_resolution_um_per_px': 0.345,  # Set to None if unknown
    'image_folder_path': "Data/images/",
}

MODEL_CONFIG = {
    'model_arc': 'Unet',
    'encoder': 'resnet34',
    'encoder_weights': 'imagenet',
    'model_folder_path': "Data/segmentation_models/",
}

TILES_SUBDIR = "tiles"

def _asymmetric_gaussian_kernel(size, sigma_left, sigma_right=None):
    """
    Generates an asymmetric Gaussian kernel for convolution.

    If sigma_left and sigma_right are not provided, it generates a standard
    symmetric Gaussian kernel.

    Args:
        size (int): The total size (number of points) of the kernel. Must be odd.
        sigma_left (float): The standard deviation for the left tail.
        sigma_right (float, optional): The standard deviation for the right tail. If None, the kernel is a symmetric gaussian with sigma = sigma_left.

    Returns:
        np.ndarray: The normalized 1D convolution kernel.
    """
    if size % 2 == 0:
        raise ValueError("Kernel size must be odd.")

    center = size // 2
    x = np.arange(size)

    if sigma_right is None:
        sigma_right = sigma_left
    
    kernel = np.zeros(size)
    
    kernel[:center] = np.exp(-(x[:center] - center)**2 / (2 * sigma_left**2))
    kernel[center] = 1.0
    kernel[center+1:] = np.exp(-(x[center+1:] - center)**2 / (2 * sigma_right**2))
    
    return kernel / np.sum(kernel)

def smear_spectrum(counts, size, sigma_left, sigma_right=None):
    """
    Applies asymmetric gaussian smearing to the track length distribution.

    Args:
        counts (np.ndarray): Track counts in the bins.
        size (int): The total size (number of points) of the kernel. Must be odd.
        sigma_left (float): The standard deviation for the left tail.
        sigma_right (float, optional): The standard deviation for the right tail. If None, the kernel is a symmetric gaussian with sigma = sigma_left.

    Returns:
        np.ndarray: The smeared track counts.
    """        
    smeared_counts = np.convolve(counts, _asymmetric_gaussian_kernel(size, sigma_left, sigma_right), mode='same')

    return smeared_counts
    

def slice_tif_to_png_tiles(images_per_group, input_dir, output_spec=TILES_SUBDIR, tile_size=IMAGE_CONFIG['img_height']):
    """
    Slices all TIFF images in a directory into PNG tiles of a specified size,
    grouping the images for sequential frame numbering.
    """
    output_dir = os.path.join(input_dir, output_spec)
    os.makedirs(output_dir, exist_ok=True)
    
    tif_files = sorted([f for f in os.listdir(input_dir) if f.lower().endswith(('.tif', '.tiff'))])
    
    print(f"Found {len(tif_files)} TIFF images to process with a group size of {images_per_group}.")

    if len(tif_files)% images_per_group != 0:
        raise ValueError(f"The total number of images ({len(tif_files)}) is not a multiple of images_per_group ({images_per_group}).")
    
    file_count = 0
    group_number = 0
    
    for filename in tif_files:
        if file_count % images_per_group == 0:
            group_number += 1
            
        input_filepath = os.path.join(input_dir, filename)
        base_name = os.path.splitext(filename)[0]
        
        try:
            img = Image.open(input_filepath)
            
            width, height = img.size
            
            x_tiles = width // tile_size
            y_tiles = height // tile_size
            
            if width % tile_size != 0:
                x_tiles += 1
            if height % tile_size != 0:
                y_tiles += 1

            print(f"Processing {filename} (Group {group_number}): {x_tiles}x{y_tiles} tiles.")
            
            tile_count = 0
            for i in range(x_tiles):
                for j in range(y_tiles):
                    left = i * tile_size
                    top = j * tile_size
                    right = min(left + tile_size, width)
                    bottom = min(top + tile_size, height)

                    box = (left, top, right, bottom)
                    
                    tile = img.crop(box)
                    
                    if tile.width == tile_size and tile.height == tile_size:
                        tile_filename = f"Frame{group_number}_{base_name}_{i}_{j}.png"
                        output_filepath = os.path.join(output_dir, tile_filename)
                    
                        tile.save(output_filepath, "PNG")
                        tile_count += 1
                    
            print(f"-> Created {tile_count} full-size tiles for {filename}")
            
            file_count += 1

        except Exception as e:
            print(f"Error processing {filename}: {e}")
            
    print(f"\n--- Slicing Complete. Total images processed: {file_count} in {group_number} groups. ---")

def convert_xml_to_masks(xml_path, output_dir):
    """
    Parses a CVAT-style XML annotation file containing multiple images and
    converts the annotations for each image into a separate binary mask file.

    Args:
        xml_path (str): The path to the input annotations.xml file.
        output_dir (str): The directory where the mask PNG files will be saved.
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Masks will be saved to: {output_dir}")

        tree = ET.parse(xml_path)
        root = tree.getroot()

        image_tags = root.findall('image')
        if not image_tags:
            print("Error: No <image> tags found in the XML file.")
            return

        for image_tag in image_tags:
            image_name = image_tag.get('name')
            img_width = int(image_tag.get('width'))
            img_height = int(image_tag.get('height'))

            final_mask = np.zeros((img_height, img_width), dtype=np.uint8)

            mask_tags = image_tag.findall('mask')
            
            for mask_tag in mask_tags:
                if mask_tag.get('label') != 'Track':
                    continue

                rle_string = mask_tag.get('rle')
                rle_parts = [int(p) for p in rle_string.split(', ')]
                
                values = []
                current_val = 0
                for run_length in rle_parts:
                    values.extend([current_val] * run_length)
                    current_val = 1 - current_val

                mask_height = int(mask_tag.get('height'))
                mask_width = int(mask_tag.get('width'))
                mask_patch = np.array(values, dtype=np.uint8).reshape(mask_height, mask_width)

                left = int(mask_tag.get('left'))
                top = int(mask_tag.get('top'))

                roi = final_mask[top : top + mask_height, left : left + mask_width]
                
                np.maximum(roi, mask_patch * 255, out=roi)
            
            base_name = os.path.splitext(image_name)[0]
            output_mask_name = f"{base_name}_mask.png"
            output_mask_path = os.path.join(output_dir, output_mask_name)

            cv2.imwrite(output_mask_path, final_mask)
          
        print(f"\nProcessed {len(image_tags)} images in total.")

    except FileNotFoundError:
        print(f"Error: The file '{xml_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def get_preprocessing(preprocessing_fn):
    """Combines model-specific normalization with tensor conversion."""
    _transform = [
        A.Lambda(image=preprocessing_fn),
        ToTensorV2(),
    ]
    return A.Compose(_transform)

class OptimusPrimus:
    def __init__(self, image_spec, model_spec, image_config=IMAGE_CONFIG, tiles_subdir=TILES_SUBDIR, model_config=MODEL_CONFIG):
        self.image_folder_path = image_config['image_folder_path']

        self.image_spec = image_spec

        self.sub_dir = tiles_subdir

        self.image_path = os.path.join(self.image_folder_path, self.image_spec, self.sub_dir)
        self.model_folder_path = model_config['model_folder_path']

        self.model_spec = model_spec
        self.model_ext = ".pth"

        self.model_path = os.path.join(self.model_folder_path, self.model_spec + self.model_ext)

        self.output_folder_path = "Data/inference_results"
        self.output_spec = image_spec + "_" + model_spec + ".csv"

        self.output_path = os.path.join(self.output_folder_path, self.output_spec)

        self.model_arc = model_config['model_arc']
        self.encoder = model_config['encoder']
        self.encoder_weights = model_config['encoder_weights']

        self.img_height = image_config['img_height']
        self.img_width = image_config['img_width']
        self.input_channels_config = image_config['input_channels_config']

        # Set to None if unknown.
        self.pixel_resolution_um_per_px = image_config['pixel_resolution_um_per_px']

        self.device = "cuda" if torch.cuda.is_available() else "cpu"

        if not os.path.isdir(self.image_path):
            raise Warning(f"ERROR: Image folder not found at '{self.image_path}'")
        if not os.path.exists(self.model_path):
            raise Warning(f"ERROR: Model file not found at '{self.model_path}'")

        print(f"Configuration Loaded:")
        print(f"  Image Folder: {self.image_path}")
        print(f"  Model Path: {self.model_path}")
        print(f"  Output Directory: {self.output_path}")
        print(f"  Using Device: {self.device}")

    def load_image_groups(self):
        if not os.path.isdir(self.image_path):
            raise ValueError(f"ERROR: Image folder not found at '{self.image_path}'")

        image_files = sorted(glob.glob(os.path.join(self.image_path, '*.png')))

        file_pattern = re.compile(r'Frame(\d+)_Acquisition_(\d+)_(\d+)_(\d+).png')

        image_groups = defaultdict(list)

        for fpath in image_files:
            fname = os.path.basename(fpath)
            match = file_pattern.match(fname)
            if match:
                n, image_id, top, left = match.groups()
                group_key = f"Frame{n}_{top}_{left}"
                image_groups[group_key].append(fpath)

        print(f"Found {len(image_files)} total images.")
        print(f"Grouped into {len(image_groups)} unique locations (Z-stacks). Each should have {len(list(image_groups.values())[0])} images.")
        img_area_cm2 = self.img_height*self.img_width*(self.pixel_resolution_um_per_px*1e-4)**2
        self.tot_area_cm2 = len(image_groups)*img_area_cm2

        return image_groups

    def load_model(self):
        try:
            model = smp.create_model(
                arch=self.model_arc,
                encoder_name=self.encoder,
                encoder_weights=None,
                in_channels=self.input_channels_config,
                classes=1,
                activation=None, 
            )
            model.load_state_dict(torch.load(self.model_path, map_location=self.device))
            model.to(self.device)
            model.eval()
            print("Model loaded successfully and set to evaluation mode.")
            return model
        except Exception as e:
            print(f"Error creating or loading model: {e}")
            raise e
        
    def perform_inference(self):

        image_groups = self.load_image_groups()
        inference_model = self.load_model()

        ellipse_analysis_results = []

        preprocessing_fn = smp.encoders.get_preprocessing_fn(self.encoder, self.encoder_weights)
        preprocessing = get_preprocessing(preprocessing_fn)

        all_keys = list(image_groups.keys())
        num_samples = 5
        selected_keys = random.sample(all_keys, min(len(all_keys), num_samples))

        for group_key, file_list in tqdm(image_groups.items(), desc="Ellipse Analysis"):
            
            group_masks = []
            original_vis_image = None
            for idx, image_path in enumerate(file_list[:]):
                image = cv2.imread(image_path)
                image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
                
                if idx == 6:
                    original_vis_image = image.copy()

                sample = preprocessing(image=image)
                preprocessed_image = sample['image'].unsqueeze(0).to(self.device, dtype=torch.float32)
                
                with torch.no_grad():
                    pred_output = inference_model(preprocessed_image)
                    pred_probs = torch.sigmoid(pred_output)
                    pr_mask = (pred_probs.squeeze().cpu().numpy().round())
                
                group_masks.append(pr_mask)
            
            aggregated_mask = np.max(np.stack(group_masks, axis=0), axis=0).astype(np.uint8)
            
            ellipse_vis_image = original_vis_image.copy()

            contours, _ = cv2.findContours(aggregated_mask.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            
            track_id_counter = 1
            
            for contour in contours:
                if cv2.contourArea(contour) < 10:
                    continue
            
                if len(contour) < 5:
                    continue
                    
                ellipse = cv2.fitEllipse(contour)

                cv2.ellipse(ellipse_vis_image, ellipse, (0, 255, 0), 2)
                
                ((center_x, center_y), (axis_a_px, axis_b_px), orientation_deg) = ellipse

                major_axis_len_px = max(axis_a_px, axis_b_px)
                minor_axis_len_px = min(axis_a_px, axis_b_px)
                
                track_props = {
                    "image_filename": os.path.basename(file_list[0]), 
                    "track_id": track_id_counter,
                    "centroid_x_px_ellipse": round(center_x, 1), 
                    "centroid_y_px_ellipse": round(center_y, 1),
                    "major_axis_px": round(major_axis_len_px, 2),
                    "minor_axis_px": round(minor_axis_len_px, 2),
                    "orientation_deg": round(orientation_deg, 2),
                }
                
                if self.pixel_resolution_um_per_px is not None:
                    track_props['len_um'] = round(major_axis_len_px * self.pixel_resolution_um_per_px, 2)

                ellipse_analysis_results.append(track_props)
                track_id_counter += 1

            if group_key in selected_keys:
                fig, axes = plt.subplots(2, 1, figsize=(6, 13))
                    
                axes[0].imshow(ellipse_vis_image)
                axes[0].set_title(f"Ellipse Fitting Overlay, example frame {group_key}")
                axes[0].axis('off')

                axes[1].imshow(aggregated_mask, cmap='gray')
                axes[1].set_title("Inferred Mask")
                axes[1].axis('off')

                plt.tight_layout()
                plt.show()

        df_ellipse_results = pd.DataFrame(ellipse_analysis_results)
        df_ellipse_results.to_csv(self.output_path, index=False)

        self.ellipses = df_ellipse_results

        return
    
    def inference_from_file(self, csv_path=None):
        if csv_path is None:
            csv_path = self.output_path

        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"ERROR: CSV file not found at '{csv_path}'")

        df_ellipse_results = pd.read_csv(csv_path)
        self.ellipses = df_ellipse_results

    def get_track_distributions(self, metric='len_um'):
        if not hasattr(self, 'ellipses'):
            raise ValueError("No ellipse data found. Please run perform_inference() first.")
        
        if not hasattr(self, 'tot_area_cm2'):
            img_area_cm2 = self.img_height*self.img_width*(self.pixel_resolution_um_per_px*1e-4)**2
            self.tot_area_cm2 = len(np.unique(self.ellipses['image_filename']))*img_area_cm2
        
        track_density = len(self.ellipses)/self.tot_area_cm2

        mean_value = self.ellipses[metric].mean()
        std_value = self.ellipses[metric].std()

        low, median, high = self.ellipses[metric].quantile([0.25, 0.5, 0.75])

        distribution_summary = {
            "track_density_per_cm2": track_density,
            "mean": mean_value,
            "std": std_value,
            "low": low,
            "median": median,
            "high": high
        }

        return distribution_summary
    
    def get_training_pairs(self, training_folder, mask_folder='training_masks', mask_extension='_mask.png'):
        image_files = []
        image_files.extend(sorted(glob.glob(os.path.join(training_folder, '*.png'))))

        file_pairs = []
        skipped_images_count = 0

        print("Attempting to pair images with masks...")
        for img_path in image_files:
            try:
                img_filename = os.path.basename(img_path)
                img_stem = os.path.splitext(img_filename)[0]

                mask_filename = img_stem + mask_extension
                mask_path = os.path.join(training_folder, mask_folder, mask_filename)

                if os.path.exists(mask_path):
                    file_pairs.append((img_path, mask_path))
                else:
                    skipped_images_count += 1

            except Exception as e:
                print(f"Error processing image path '{img_path}': {e}")
                skipped_images_count += 1

        print(f"\nSuccessfully created {len(file_pairs)} image-mask pairs.")
        if skipped_images_count > 0:
            print(f"Warning: Skipped {skipped_images_count} images because their corresponding masks were not found.")

        return file_pairs

    

    def evaluate_binned_efficiency(self, training_folder, mask_folder='training_masks', mask_extension='_mask.png', iou_threshold=0.5, num_bins=10):
        """
        Calculates size-binned Precision and Recall by matching predicted 
        instances to ground truth instances via IoU.
        """
        validation_pairs = self.get_training_pairs(training_folder, mask_folder, mask_extension)
        inference_model = self.load_model()
        
        preprocessing_fn = smp.encoders.get_preprocessing_fn(self.encoder, self.encoder_weights)
        preprocessing = get_preprocessing(preprocessing_fn)

        gt_log = []
        pred_log = []

        for img_path, mask_path in tqdm(validation_pairs, desc="Evaluating Efficiency"):
            image = cv2.imread(img_path)
            image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
            gt_mask_full = cv2.imread(mask_path, cv2.IMREAD_GRAYSCALE)
            
            sample = preprocessing(image=image_rgb)
            preprocessed_image = sample['image'].unsqueeze(0).to(self.device, dtype=torch.float32)
            with torch.no_grad():
                pred_output = inference_model(preprocessed_image)
                pred_probs = torch.sigmoid(pred_output)
                pred_mask_full = (pred_probs.squeeze().cpu().numpy().round()).astype(np.uint8)

            def extract_instances(binary_mask):
                contours, _ = cv2.findContours(binary_mask.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                instances = []
                for idx, contour in enumerate(contours):
                    if cv2.contourArea(contour) < 10 or len(contour) < 5:
                        continue
                    
                    ellipse = cv2.fitEllipse(contour)
                    ((_, _), (axis_a_px, axis_b_px), _) = ellipse
                    major_axis_len_px = max(axis_a_px, axis_b_px)
                    
                    size_metric = major_axis_len_px
                    if self.pixel_resolution_um_per_px is not None:
                        size_metric = major_axis_len_px * self.pixel_resolution_um_per_px

                    instance_mask = np.zeros_like(binary_mask)
                    cv2.drawContours(instance_mask, [contour], -1, 1, thickness=cv2.FILLED)
                    
                    instances.append({'id': idx, 'size': size_metric, 'mask': instance_mask, 'matched': False})
                return instances

            gt_instances = extract_instances(gt_mask_full)
            pred_instances = extract_instances(pred_mask_full)

            for pred in pred_instances:
                best_iou = 0
                best_gt_idx = -1
                
                for idx, gt in enumerate(gt_instances):
                    if gt['matched']:
                        continue
                        
                    intersection = np.logical_and(pred['mask'], gt['mask']).sum()
                    if intersection == 0: continue
                    
                    union = np.logical_or(pred['mask'], gt['mask']).sum()
                    iou = intersection / union
                    
                    if iou > best_iou:
                        best_iou = iou
                        best_gt_idx = idx
                        
                if best_iou >= iou_threshold:
                    pred['matched'] = True
                    gt_instances[best_gt_idx]['matched'] = True

            for gt in gt_instances:
                gt_log.append({'len_um': gt['size'], 'is_true_positive': gt['matched']})
            for pred in pred_instances:
                pred_log.append({'len_um': pred['size'], 'is_true_positive': pred['matched']})

        df_gt = pd.DataFrame(gt_log)
        df_pred = pd.DataFrame(pred_log)

        self.gt_ellipses = df_gt
        self.pred_ellipses = df_pred

        if df_gt.empty or df_pred.empty:
            print("Warning: No instances found in Ground Truth or Predictions.")
            return None

        df_gt['size_bin'], bins = pd.qcut(df_gt['len_um'], q=num_bins, retbins=True, duplicates='drop')
        
        df_pred['size_bin'] = pd.cut(df_pred['len_um'], bins=bins, include_lowest=True)

        recall_df = df_gt.groupby('size_bin')['is_true_positive'].agg(
            Recall='mean',
            GT_Count='count'
        )

        precision_df = df_pred.groupby('size_bin')['is_true_positive'].agg(
            Precision='mean', 
            Pred_Count='count'
        )

        bins_mid = (bins[:-1] + bins[1:]) / 2

        efficiency_table = pd.concat([recall_df, precision_df], axis=1)
        efficiency_table.index.name = f'Size Bin ({ "um" if self.pixel_resolution_um_per_px else "px" })'

        efficiency_table['Bin_mid'] = bins_mid

        efficiency_table = efficiency_table.fillna(0)

        efficiency_table.to_csv(os.path.join(self.output_folder_path, f"{self.model_spec}_efficiency_by_size_bin.csv"))

        self.efficiency_table = efficiency_table

        return efficiency_table
    
    def efficiency_distribution_from_file(self, csv_path=None):
        if csv_path is None:
            csv_path = os.path.join(self.output_folder_path, f"{self.model_spec}_efficiency_by_size_bin.csv")

        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"ERROR: CSV file not found at '{csv_path}'")

        efficiency_table = pd.read_csv(csv_path, index_col=0)
        self.efficiency_table = efficiency_table

    def get_training_set_distributions(self, metric='len_um'):
        if not hasattr(self, 'gt_ellipses'):
            raise ValueError("No training ellipse data found. Please run training_set_analysis() first.")
        
        img_area_cm2 = self.img_height*self.img_width*(self.pixel_resolution_um_per_px*1e-4)**2
        self.training_tot_area_cm2 = len(np.unique(self.gt_ellipses['image_filename']))*img_area_cm2
        
        training_track_density = len(self.gt_ellipses)/self.training_tot_area_cm2

        mean_value = self.gt_ellipses[metric].mean()
        std_value = self.gt_ellipses[metric].std()

        low, median, high = self.gt_ellipses[metric].quantile([0.25, 0.5, 0.75])

        distribution_summary = {
            "track_density_per_cm2": training_track_density,
            "mean": mean_value,
            "std": std_value,
            "low": low,
            "median": median,
            "high": high
        }

        return distribution_summary


    def apply_detection_model_efficiency(self, x_bins, counts, meas_error=1000.):
        """
        Multiplies counts (sliced) with the efficiency function for the detection model
        
        Args:
            x_bins (np.ndarray): The bin edges for the track length spectrum R [nm].
            counts (np.ndarray): The array of track counts N(R) in each bin.

        Returns:
            np.ndarray: Efficiency-corrected counts.
        """

        x_mids = x_bins[:-1] + np.diff(x_bins) / 2.0

        recall = np.interp(x_mids, self.efficiency_table['Bin_mid'], self.efficiency_table['Recall'])
        precision = np.interp(x_mids, self.efficiency_table['Bin_mid'], self.efficiency_table['Precision'])

        counts_with_efficiency = counts * recall / precision

        counts_with_measure = smear_spectrum(counts_with_efficiency, len(x_bins)//2*2-1, meas_error/np.diff(x_bins)[0], meas_error/np.diff(x_bins)[0])

        return counts_with_measure
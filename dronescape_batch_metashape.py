#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch script to process multiple imagery directories with dronescape_metashape.py

This script processes multiple sites/dates by running the main processing script
for each imagery directory in the list below.

Author: Juan C. Montes Herrera (University of Tasmania)
"""

import sys
from pathlib import Path
import Metashape

# Import the main processing function
from dronescape_metashape import main

# =============================================================================
# CONFIGURATION - Edit these variables for your batch processing
# =============================================================================

# Output directory where all Metashape projects will be saved
OUTPUT_DIR = r"F:\metashapes"

# List of imagery directories to process
# Each path should point to a YYYYMMDD/imagery/ directory

IMAGERY_DIRS = [
    # r"D:\TERN-Dronescape\NTABRT0001\20240829\imagery",
    # r"D:\TERN-Dronescape\NTABRT0002\20240829\imagery",
    # r"D:\TERN-Dronescape\NTABRT0003\20240829\imagery",
    # r"D:\TERN-Dronescape\NTABRT0004\20240828\imagery",
    # r"D:\TERN-Dronescape\NTABRT0005\20240828\imagery",
    # r"D:\TERN-Dronescape\NTABRT0006\20240828\imagery",
    # r"D:\TERN-Dronescape\NTABRT0007\20240828\imagery",
    # r"D:\TERN-Dronescape\NTABRT0009\20240829\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0001\20240902\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0002\20240901\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0003\20240902\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0004\20240902\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0005\20240902\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0006\20240902\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0007\20240831\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0008\20240831\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0009\20240901\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0010\20240901\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0011\20240901\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0012\20240924\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0013\20240925\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0014\20240925\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0017\20240924\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0018\20240924\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0019\20240924\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0024\20240926\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0025\20240926\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0026\20240830\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0027\20240830\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0028\20240927\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0029\20240927\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0030\20240927\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0031\20240927\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0032\20240830\imagery",
    # r"D:\TERN-Dronescape\NTAFIN0033\20240830\imagery",
    # r"D:\TERN-Dronescape\NTAMAC0001\20240831\imagery",
    # r"D:\TERN-Dronescape\NTAMAC0002\20240831\imagery",
    # r"D:\TERN-Dronescape\NTAMAC0003\20240831\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0004\20240722\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0005\20240722\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0006\20240723\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0007\20240721\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0008\20240723\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0009\20240720\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0010\20240721\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0011\20240723\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0012\20240722\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0013\20240718\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0014\20240720\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0015\20240721\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0017\20240720\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0018\20240723\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0019\20240722\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0020\20240717\imagery",
    # r"D:\TERN-Dronescape\SAAEYB0037\20240717\imagery",
    # r"D:\TERN-Dronescape\SAAGAW0004\20241001\imagery",
    # r"D:\TERN-Dronescape\SAAGAW0005\20241001\imagery",
    # r"D:\TERN-Dronescape\SAAGAW0006\20241002\imagery",
    # r"D:\TERN-Dronescape\SAAGAW0007\20241002\imagery",
    # r"D:\TERN-Dronescape\SAAGAW0008\20241001\imagery",
    # r"D:\TERN-Dronescape\SAAGAW0009\20241001\imagery",
    # r"D:\TERN-Dronescape\SAAGVD0005\20240928\imagery",
    # r"D:\TERN-Dronescape\SAAGVD0006\20240928\imagery",
    # r"D:\TERN-Dronescape\SAAKAN0001\20241208\imagery",
    # r"D:\TERN-Dronescape\SAAKAN0002\20241210\imagery",
    # r"D:\TERN-Dronescape\SAAKAN0003\20241209\imagery",
    # r"D:\TERN-Dronescape\SAAKAN0004\20241210\imagery",
    # r"D:\TERN-Dronescape\SAAKAN0005\20241210\imagery",
    # r"D:\TERN-Dronescape\SAAKAN0006\20241211\imagery",
    # r"D:\TERN-Dronescape\SAAKAN0007\20241211\imagery",
    # r"D:\TERN-Dronescape\SAAKAN0008\20241208\imagery",
    # r"D:\TERN-Dronescape\SAAKAN0009\20241208\imagery",
    # r"D:\TERN-Dronescape\SAAKAN0010\20241207\imagery",
    # r"D:\TERN-Dronescape\SAAKAN0011\20241207\imagery",
    # r"D:\TERN-Dronescape\SAAKAN0012\20241207\imagery",
    # r"D:\TERN-Dronescape\SAAKAN0013\20241211\imagery",
    # r"D:\TERN-Dronescape\SAANCP0001\20241122\imagery",
    # r"D:\TERN-Dronescape\SAANCP0002\20241121\imagery",
    # r"D:\TERN-Dronescape\SAANCP0003\20241121\imagery",
    # r"D:\TERN-Dronescape\SAANCP0004\20241121\imagery",
    # r"D:\TERN-Dronescape\SAANCP0005\20241121\imagery",
    # r"D:\TERN-Dronescape\SAANCP0006\20241120\imagery",
    # r"D:\TERN-Dronescape\SAANCP0007\20241120\imagery",
    # r"D:\TERN-Dronescape\SAANCP0008\20241118\imagery",
    # r"D:\TERN-Dronescape\SAANCP0009\20241119\imagery",
    # r"D:\TERN-Dronescape\SAASTP0033\20241002\imagery",
    # r"D:\TERN-Dronescape\SAASTP0034\20241002\imagery",
    # r"D:\TERN-Dronescape\SAASTP0036\20240928\imagery",
    # r"D:\TERN-Dronescape\SAASTP0037\20240929\imagery",
    # r"D:\TERN-Dronescape\SAASVP0001\20241120\imagery",
    # r"D:\TERN-Dronescape\SAASVP0002\20241120\imagery",
    # r"D:\TERN-Dronescape\SASMDD0001\20240312\imagery",
    # r"D:\TERN-Dronescape\SASMDD0002\20240312\imagery",
    # r"D:\TERN-Dronescape\SASMDD0003\20240312\imagery",
    # r"D:\TERN-Dronescape\SASMDD0008\20240313\imagery",
    # r"D:\TERN-Dronescape\SASMDD0009\20240313\imagery",
    # r"D:\TERN-Dronescape\SASMDD0010_18\20240311\imagery",
    # r"D:\TERN-Dronescape\SASMDD0011\20240314\imagery",
    # r"D:\TERN-Dronescape\SASMDD0012\20240314\imagery",
    # r"D:\TERN-Dronescape\SASMDD0013\20240310\imagery",
    # r"D:\TERN-Dronescape\SASMDD0014\20240313\imagery",
    # r"D:\TERN-Dronescape\SASMDD0016_17\20240315\imagery",
    # r"D:\TERN-Dronescape\SASMDD005_6\20240314\imagery",
    # r"D:\TERN-Dronescape\SASRIV0001\20240313\imagery",
    # r"D:\TERN-Dronescape\SASRIV0002\20240315\imagery",
    # r"D:\TERN-Dronescape\WAAAVW0003\20250405\imagery",
    # r"D:\TERN-Dronescape\WAAAVW0004\20250405\imagery",
    # r"D:\TERN-Dronescape\WAAAVW0008\20250405\imagery",
    # r"D:\TERN-Dronescape\WAACOO0020\20250823\imagery",
    # r"D:\TERN-Dronescape\WAACOO0021\20250822\imagery",
    # r"D:\TERN-Dronescape\WAACOO0022\20250822\imagery",
    # r"D:\TERN-Dronescape\WAACOO0023\20250822\imagery",
    # r"D:\TERN-Dronescape\WAACOO0024\20250824\imagery",
    # r"D:\TERN-Dronescape\WAACOO0025\20250824\imagery",
    r"D:\TERN-Dronescape\WAACOO0034\20250827\imagery",
    # r"D:\TERN-Dronescape\WAACOO0039\20250823\imagery",
    # r"D:\TERN-Dronescape\WAACOO0040\20250823\imagery",
    # r"D:\TERN-Dronescape\WAACOO0041\20250823\imagery",
    # r"D:\TERN-Dronescape\WAACOO0042\20250822\imagery",
    # r"D:\TERN-Dronescape\WAACOO0043\20250822\imagery",
    # r"D:\TERN-Dronescape\WAACOO0044\20250822\imagery",
    # r"D:\TERN-Dronescape\WAAESP0001\20250328\imagery",
    # r"D:\TERN-Dronescape\WAAJAF0001\20250329\imagery",
    # r"D:\TERN-Dronescape\WAAJAF0002\20250329\imagery",
    # r"D:\TERN-Dronescape\WAAJAF0003\20250328\imagery",
    # r"D:\TERN-Dronescape\WAAJAF0007\20250404\imagery",
    # r"D:\TERN-Dronescape\WAAJAF0009\20250404\imagery",
    # r"D:\TERN-Dronescape\WAAJAF0010\20250404\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0001\20250529\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0002\20250701\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0003\20250701\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0004\20250702\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0005\20250701\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0006\20250701\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0007\20250701\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0008\20250630\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0009\20250630\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0010\20250603\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0011\20250603\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0012\20250529\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0014\20250529\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0015\20250601\imagery",
    r"D:\TERN-Dronescape\WAAPIL0016\20250530\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0017\20250629\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0018\20250528\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0019\20250528\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0020\20250702\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0021\20250702\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0022\20250702\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0023\20250629\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0024\20250629\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0025\20250726\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0026\20250726\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0027\20250726\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0028\20250726\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0029\20250726\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0030\20250727\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0031\20250728\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0032\20250728\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0033\20250728\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0034\20250727\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0036\20250602\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0037\20250601\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0038\20250601\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0039\20250602\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0040\20250531\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0041\20250531\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0042\20250530\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0042\20250531\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0043\20250604\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0044\20250604\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0045\20250604\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0046\20250723\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0047\20250723\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0048\20250724\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0049\20250723\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0050\20250723\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0051\20250723\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0052\20250625\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0053\20250625\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0054\20250626\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0055\20250626\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0056\20250626\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0057\20250625\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0058\20250625\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0059\20250626\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0060\20250627\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0061\20250627\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0062\20250627\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0063\20250627\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0064\20250628\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0065\20250628\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0066\20250628\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0071\20250725\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0072\20250725\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0073\20250725\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0074\20250725\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0075\20250725\imagery",
    # r"D:\TERN-Dronescape\WAAPIL0076\20250725\imagery",
    # r"D:\TERN-Dronescape\WAASWA0012\20250326\imagery",
    # r"D:\TERN-Dronescape\WAASWA0013\20250326\imagery",
    # r"D:\TERN-Dronescape\WAASWA0014\20250326\imagery",
    # r"D:\TERN-Dronescape\WAFWAR0005\20250329\imagery",
    # r"D:\TERN-Dronescape\WAFWAR0007\20250330\imagery",
    # r"D:\TERN-Dronescape\WAFWAR0008\20250330\imagery",
    # r"D:\TERN-Dronescape\WAGCOO0004\20250821\imagery"
]

# Optional: Enable oblique cameras - set to True to enable
ENABLE_OBLIQUE = False

# =============================================================================
# BATCH PROCESSING - No need to edit below this line
# =============================================================================

def batch_process():
    """Process all imagery directories in the list."""
    
    total = len(IMAGERY_DIRS)
    print(f"Starting batch processing of {total} directories")
    print("=" * 80)
    
    for i, imagery_dir in enumerate(IMAGERY_DIRS, 1):
        print(f"\n[{i}/{total}] Processing: {imagery_dir}")
        print("-" * 80)
        
        # Verify directory exists
        if not Path(imagery_dir).exists():
            print(f"WARNING: Directory does not exist, skipping: {imagery_dir}")
            continue
        
        # Build command line arguments
        args = [
            "dronescape_metashape.py",
            "-imagery_dir", imagery_dir,
            "-out", OUTPUT_DIR,
        ]
        
        # Add oblique flag if enabled
        if ENABLE_OBLIQUE:
            args.append("-enable_oblique")
        
        # Temporarily replace sys.argv to simulate command line arguments
        original_argv = sys.argv
        sys.argv = args
        
        try:
            # Run the main processing function
            main()
            print(f"\n[{i}/{total}] Successfully completed: {imagery_dir}")
        except Exception as e:
            print(f"\n[{i}/{total}] ERROR processing {imagery_dir}: {e}")
            print("Continuing to next directory...")
            
            # Cleanup on error to ensure next run starts fresh
            try:
                doc = Metashape.app.document
                doc.clear()
                Metashape.app.releaseFreeMemory()
                print("Performed emergency cleanup after error.")
            except Exception as cleanup_error:
                print(f"Warning: Cleanup after error failed: {cleanup_error}")
        finally:
            # Restore original sys.argv
            sys.argv = original_argv
            
            # Release free memory between runs
            try:
                Metashape.app.releaseFreeMemory()
            except Exception as cleanup_error:
                print(f"Warning: Memory cleanup failed: {cleanup_error}")
        
        print("=" * 80)
    
    print(f"\nBatch processing complete! Processed {total} directories.")


if __name__ == "__main__":
    batch_process()


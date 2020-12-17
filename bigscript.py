import pandas as pd
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import photutils
from photutils import segmentation as seg
import scipy.ndimage as ndi
import statmorph
from astropy.visualization import (AsinhStretch, LogStretch, ImageNormalize)
import sys
from astropy import wcs
import urllib.request
from urllib.error import HTTPError, URLError
import multiprocessing as mp
import os
from os import path
import time
import statistics
from matplotlib import gridspec
import math


# Delete the temporary galaxy files
def clearTempFiles(galaxy):
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G.fits'):
        os.system(f'rm /mnt/scratch/temp_galaxy_storage/{galaxy}_G.fits')  
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_model.fits'):
        os.system(f'rm /mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_model.fits')
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_model.fits'):
        os.system(f'rm /mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_model.fits') 
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_mask.fits'):
        os.system(f'rm /mnt/scratch/temp_galaxy_storage/{galaxy}_G_mask.fits')
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_psf.fits'):
        os.system(f'rm /mnt/scratch/temp_galaxy_storage/{galaxy}_G_psf.fits')
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_sig.fits'):
        os.system(f'rm /mnt/scratch/temp_galaxy_storage/{galaxy}_G_sig.fits')
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_residual.fits'):
        os.system(f'rm /mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_residual.fits')
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_residual.fits'):
        os.system(f'rm /mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_residual.fits')


# Check if any galaxy files downloaded are 0 MB
def checkCorrupt(galaxy):
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G.fits'):
        if(os.path.getsize(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G.fits') == 0):
            return True;
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_model.fits'):
        if(os.path.getsize(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_model.fits') == 0):
            return True;
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_model.fits'):
        if(os.path.getsize(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_model.fits') == 0):
            return True;
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_mask.fits'):
        if(os.path.getsize(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_mask.fits') == 0):
            return True;
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_psf.fits'):
        if(os.path.getsize(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_psf.fits') == 0):
            return True;
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_sig.fits'):
        if(os.path.getsize(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_sig.fits') == 0):
            return True;
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_residual.fits'):
        if(os.path.getsize(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_residual.fits') == 0):
            return True;
    if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_residual.fits'):
        if(os.path.getsize(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_residual.fits') == 0):
            return True;
    return False;

    
# Determine if pixel [i,j] is within ellipse area
def inEllipse(i,j,xCenter,yCenter,size,ratio,angle):

    # Size is effective radius of galaxy which is then multiplied by a number to increase size of ellipse
    if (size > 30):
        size = 20*size
    else:   
        size = 10*size

    x = (i-xCenter)
    y = (j-yCenter)
    a = math.sqrt((size**2)/ratio)
    b = (size**2)/a

    if(  (((((x*math.cos(math.radians(angle))) - (y*math.sin(math.radians(angle))))**2)  / (a**2)) + ((((y*math.cos(math.radians(angle))) + (x*math.sin(math.radians(angle))))**2)  / (b**2))) <= 1):
        return True


def statmorphWrapper(index_pairs):

    crash = ['NGVSJ12:29:48.87+13:25:46.0', 'NGVSJ12:30:49.42+12:23:28.0']
    
    df = pd.read_csv('NGVSgalaxies.csv')

    # Iterate through rows in csv file containing measurements for each galaxy
    for row in df.iloc[index_pairs[0]:index_pairs[1]].itertuples(index=True, name='Pandas'):
        galaxy = row.Official_name
        base = 'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/files/vault/ngvs/data/NGVS/galaxies/'
        galaxyPath = f'{galaxy}/{galaxy}_G'

        if galaxy in crash :
            continue

        startcopy = time.time()
        # Try to copy galaxy files
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxy}/{galaxy}_G.fits /mnt/scratch/temp_galaxy_storage')
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxyPath}_iso_model.fits /mnt/scratch/temp_galaxy_storage')
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxyPath}_galfit_model.fits /mnt/scratch/temp_galaxy_storage')
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxyPath}_mask.fits /mnt/scratch/temp_galaxy_storage')
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxyPath}_psf.fits /mnt/scratch/temp_galaxy_storage')
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxyPath}_sig.fits /mnt/scratch/temp_galaxy_storage')
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxyPath}_iso_residual.fits /mnt/scratch/temp_galaxy_storage')
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxyPath}_galfit_residual.fits /mnt/scratch/temp_galaxy_storage')
        endcopy = time.time() - startcopy

        # If one of the required files is missing, continue to next galaxy                                                           
        if(not path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G.fits') or (not path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_model.fits') and not \
        path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_model.fits')) or not path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_mask.fits') or not \
        path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_psf.fits') or not path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_sig.fits') or (not \
        path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_residual.fits') and not path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_residual.fits'))):
            print(f'missing {galaxy}')
            writeFile = open(f'/mnt/scratch/missing/{galaxy}.txt', 'w')     
            writeFile.write(f'{galaxy}\n')
            if(not path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G.fits')):
                writeFile.write('missing original\n')
            if(not path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_model.fits')):
                writeFile.write('missing iso model\n')
            if(not path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_model.fits')):
                writeFile.write('missing galfit model\n')
            if(not path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_mask.fits')):
                writeFile.write('missing mask\n')
            if(not path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_psf.fits')):
                writeFile.write('missing psf\n')
            if(not path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_sig.fits')):
                writeFile.write('missing sig\n')
            if(not path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_residual.fits')):
                writeFile.write('missing iso residual\n')
            if(not path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_residual.fits')):
                writeFile.write('missing galfit residual\n')
            
            writeFile.close()
            clearTempFiles(galaxy)
            continue

        # ONLY PROCESS LARGE FILE
        #---------------------------------------------------------------------------------------
        if(os.path.getsize(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G.fits') < 300000000):
            clearTempFiles(galaxy)
            continue
        #---------------------------------------------------------------------------------------

        # Create check file for current galaxy
        print(f'checking {galaxy}')
        writeFile = open(f'/mnt/scratch/check/{galaxy}.txt', 'w') 
        writeFile.write('checking')
        writeFile.close()


        # If any of the galaxy files are empty then create galaxy corrupt file and continue to next galaxy
        if(checkCorrupt(galaxy)):
            writeFile = open(f'/mnt/scratch/corrupt/{galaxy}.txt', 'w')     
            writeFile.write(f'corrupt\n')
            writeFile.close()
            clearTempFiles(galaxy)
            continue

        # Beginning of segmentation map creation

        startseg = time.time()

        hdu = fits.open(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G.fits')
        im_header = hdu[0].header
        im_data = hdu[0].data

         # Sky subtraction from original image
        sky_data = np.zeros(np.shape(im_data))
        sky_data += row.SKY
        im_sky_subtracted = im_data - sky_data

        # Calculate nucleus center
        ra = row.NGVS_ra
        dec = row.NGVS_dec
        mywcs = wcs.WCS(im_header)
        xCenter, yCenter = mywcs.all_world2pix([[ra, dec]], 0)[0]

        mask_data = fits.getdata(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_mask.fits')

        # If a iso model exists then use that file for the original model data, otherwise use galfit
        if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_model.fits'):
            original_model_data = fits.getdata(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_model.fits')
        else:
            original_model_data = fits.getdata(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_model.fits')

        # If nucleus exists then mask nucleus
        if(row.Nuc_Flag == 1):
            for i in range(len(mask_data)):
                for j in range(len(mask_data)):
                    # TODO: Change radius of nucleus mask
                    if(((i-xCenter)**2) + ((j-yCenter)**2) <= (5**2)):
                        mask_data[i][j] = 100
                    
        
        ellipse_data = np.zeros(np.shape(original_model_data))
        # Calculate median of original model values within 10 pixel radius from nucleus center
        pixelList = []
        for i in range(len(original_model_data)):
            for j in range(len(original_model_data)):
                    if(((i-xCenter)**2) + ((j-yCenter)**2) <= (10**2) and mask_data[i][j] != 100):
                        pixelList.append(original_model_data[i][j])   

        median = statistics.median(pixelList)

        # Create Segmentation Map
        seg_data = np.zeros(np.shape(original_model_data))
        ellipse_data = np.zeros(np.shape(original_model_data))
        
        # isEmpty flag is used for checking if a segmentation map is valid for processing. If segmentation map 2D list is all 0's then script crashes.
        isEmpty = True

        # if median is greater than 2*sky value then create segmentation map from original model values greater than 1.4*sky value within ellipse area
        if(median > 2*row.SKY):
            for i in range(len(original_model_data)):   
                for j in range(len(original_model_data)):
                        if(inEllipse(i,j,xCenter,yCenter,row.Size,row.AxisRatio,row.PA)):
                            ellipse_data[i][j] = 100
                            if(original_model_data[i][j] > (1.4*row.SKY)):
                                seg_data[i][j] = 100
                                isEmpty = False
        # If median is less than 2*sky value then create segmentation map from original model values greater than 1.1*sky value within ellipse area
        else:
            for i in range(len(original_model_data)):
                for j in range(len(original_model_data)):
                        if(inEllipse(i,j,xCenter,yCenter,row.Size,row.AxisRatio,row.PA)):
                            ellipse_data[i][j] = 100
                            if(original_model_data[i][j] > (1.1*row.SKY)):
                                seg_data[i][j] = 100
                                isEmpty = False

        psf = fits.getdata(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_psf.fits')
        weightmap = fits.getdata(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_sig.fits')
        mask_data = np.array(mask_data, dtype=bool)

        endseg = time.time() - startseg
 
        # End of segmentation map creation


        # If the galaxy's segmentation map is empty with no area of interest, then create empty galaxy file and continue to next galaxy
        if(isEmpty):
            writeFile = open(f'/mnt/scratch/emptyseg/{galaxy}.txt', 'w')     
            writeFile.write('empty')
            writeFile.close()
            clearTempFiles(galaxy)
            continue

        start_time = time.time()

        # run statmorph on current galaxy
        source_morphs = statmorph.source_morphology(im_sky_subtracted, seg_data, mask=mask_data, weightmap=weightmap, psf=psf)
        end_time = time.time() - start_time

        morph = source_morphs[0]

        startmodelcreate = time.time()

        # create model from statmorph results
        ny, nx = im_sky_subtracted.shape
        y, x = np.mgrid[0:ny, 0:nx]
        fitted_model = statmorph.ConvolvedSersic2D(
            amplitude=morph.sersic_amplitude,
            r_eff=morph.sersic_rhalf,
            n=morph.sersic_n,
            x_0=morph.sersic_xc,
            y_0=morph.sersic_yc,
            ellip=morph.sersic_ellip,
            theta=morph.sersic_theta)
        fitted_model.set_psf(psf)
        output_model_data = fitted_model(x, y)

        endmodelcreate = time.time() - startmodelcreate
        
        if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_residual.fits'):
            original_res_data = fits.getdata(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_residual.fits')
        else:
            original_res_data = fits.getdata(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_residual.fits')


        startfig = time.time()

        # normalize images, models, segmentation map
        output_res_data = im_sky_subtracted - output_model_data
        p1 = 10. ; p2 = 90.

        im_p1 = np.percentile(im_sky_subtracted.ravel(), p1)
        im_p2 = np.percentile(im_sky_subtracted.ravel(), p2)
        normSky = ImageNormalize(im_sky_subtracted, vmin=im_p1, vmax=im_p2)

        im_p1 = np.percentile(output_model_data.ravel(), p1)
        im_p2 = np.percentile(output_model_data.ravel(), p2)
        normOutputMod = ImageNormalize(output_model_data, vmin=im_p1, vmax=im_p2)

        im_p1 = np.percentile(original_model_data.ravel(), p1)
        im_p2 = np.percentile(original_model_data.ravel(), p2)
        normOriginalMod = ImageNormalize(original_model_data, vmin=im_p1, vmax=im_p2)

        im_p1 = np.percentile(output_res_data.ravel(), p1)
        im_p2 = np.percentile(output_res_data.ravel(), p2)
        normOutputRes = ImageNormalize(output_res_data, vmin=im_p1, vmax=im_p2)

        im_p1 = np.percentile(original_res_data.ravel(), p1)
        im_p2 = np.percentile(original_res_data.ravel(), p2)
        normOriginalRes = ImageNormalize(original_res_data, vmin=im_p1, vmax=im_p2)

        # create figures for images, models, segmentation map
        gs = gridspec.GridSpec(2, 4, width_ratios=[1, 1, 1, 1],
         wspace=0.2, hspace=0, top=0.7, bottom=0.05, left=0.1, right=0.5)

        fig = plt.figure(figsize=(30,10))

        ax= plt.subplot(gs[0,0])
        ax.imshow(im_sky_subtracted, norm=normSky, cmap='gray', origin='lower')
        ax.set_title('Sky Subtracted Image', fontsize=15)

        ax= plt.subplot(gs[0,1])
        ax.imshow(original_model_data, norm=normOriginalMod, cmap='gray', origin='lower')
        ax.set_title('Original Model', fontsize=15)

        ax= plt.subplot(gs[0,2])
        ax.imshow(output_model_data, norm=normOutputMod, cmap='gray', origin='lower')
        ax.set_title('Output Model', fontsize=15)

        ax= plt.subplot(gs[0,3])
        ax.imshow(mask_data, cmap='gray', origin='lower')
        ax.set_title('Mask', fontsize=15)

        ax= plt.subplot(gs[1,0])
        ax.imshow(seg_data, cmap='gray', origin='lower')
        ax.set_title('Segmap', fontsize=15)

        ax= plt.subplot(gs[1,1])
        ax.imshow(ellipse_data, cmap='gray', origin='lower')
        ax.set_title('Ellipse Area', fontsize=15)

        ax= plt.subplot(gs[1,2])
        ax.imshow(original_res_data, norm=normOriginalRes, cmap='gray', origin='lower')
        ax.set_title('Original Residual', fontsize=15)

        ax= plt.subplot(gs[1,3])
        ax.imshow(output_res_data, norm=normOutputRes, cmap='gray', origin='lower')
        ax.set_title('Output Residual', fontsize=15)

        endfig = time.time() - startfig

        # save figures as PNG image to output directory
        fig.savefig(f'/mnt/scratch/output/{galaxy}_sourcemorph:{round(end_time, 2)}_seg={round(endseg, 2)}_RE={round(row.Size, 3)}_mag={round(row.principleg_mag_cg, 3)}.png', facecolor='w', edgecolor='w', transparent=False, bbox_inches='tight')
        plt.close(fig)

        
        # UNCOMMENT TO SAVE AS MULTI EXTENSION FITS FILE INSTEAD OF PNG 
        #------------------------------------------------------------------------------------------------------------------------ 
        # primary_hdu = fits.PrimaryHDU(im_sky_subtracted, header=im_header)
        # image_hdu = fits.ImageHDU(output_model_data)
        # image_hdu2 = fits.ImageHDU(output_res_data)
        # hdul = fits.HDUList([primary_hdu, image_hdu, image_hdu2])

        # upload fits file to VOSpace
        # hdul.writeto(f'/mnt/scratch/output/{galaxy}_output.fits', overwrite=True)
        # os.system(f'vcp /mnt/scratch/output/{galaxy}_output.fits vos:ngvs/data/STATMORPH/FITS_output/{galaxy}_output.fits')
        # if path.isfile(f'/mnt/scratch/output/{galaxy}_output.fits'):
        #    os.system(f'rm /mnt/scratch/output/{galaxy}_output.fits')
        #------------------------------------------------------------------------------------------------------------------------
        
        # UPLOAD PNG FILE WITH FLAGS
        # os.system(f'vcp /mnt/scratch/output/{galaxy}_time:{round(end_time, 2)}_Flag={morph.flag}_SersicFlag={morph.flag_sersic}.png \
        # vos:ngvs/data/STATMORPH/filesize_bug/{galaxy}_time:{round(end_time, 2)}_Flag={morph.flag}_SersicFlag={morph.flag_sersic}.png')
        # if path.isfile(f'/mnt/scratch/output/{galaxy}_time:{round(end_time, 2)}_Flag={morph.flag}_SersicFlag={morph.flag_sersic}.png'):
        #    os.system(f'rm /mnt/scratch/output/{galaxy}_time:{round(end_time, 2)}_Flag={morph.flag}_SersicFlag={morph.flag_sersic}.png')


        # UPLOAD PNG FILE WITH MEDIAN & SKY VALUES
        # os.system(f'vcp /mnt/scratch/output/{galaxy}_time:{round(end_time, 2)}_size={row.Size}_median={median}_2*sky={2*row.SKY}sky={row.SKY}.png \
        # vos:ngvs/data/STATMORPH/memory_bug/{galaxy}_time:{round(end_time, 2)}_size={row.Size}_median={median}_2*sky={2*row.SKY}sky={row.SKY}.png')
        # if path.isfile(f'/mnt/scratch/output/{galaxy}_time:{round(end_time, 2)}_size={row.Size}_median={median}_2*sky={2*row.SKY}sky={row.SKY}.png'):
        #    os.system(f'rm /mnt/scratch/output/{galaxy}_time:{round(end_time, 2)}_size={row.Size}_median={median}_2*sky={2*row.SKY}sky={row.SKY}.png')

        # UPLOAD PNG FILE WITH RUNNING TIMES & RE FACTOR & MAGNITUDE
        os.system(f'vcp /mnt/scratch/output/{galaxy}_sourcemorph:{round(end_time, 2)}_seg={round(endseg, 2)}_RE={round(row.Size, 3)}_mag={round(row.principleg_mag_cg, 3)}.png \
        vos:ngvs/data/STATMORPH/memory_fix/{galaxy}_sourcemorph:{round(end_time, 2)}_seg={round(endseg, 2)}_RE={round(row.Size, 3)}_mag={round(row.principleg_mag_cg, 3)}.png')
        if path.isfile(f'/mnt/scratch/output/{galaxy}_sourcemorph:{round(end_time, 2)}_seg={round(endseg, 2)}_RE={round(row.Size, 3)}_mag={round(row.principleg_mag_cg, 3)}.png'):
            os.system(f'rm /mnt/scratch/output/{galaxy}_sourcemorph:{round(end_time, 2)}_seg={round(endseg, 2)}_RE={round(row.Size, 3)}_mag={round(row.principleg_mag_cg, 3)}.png')

        hdu.close()
        clearTempFiles(galaxy)
        print(f'complete {galaxy}')
        
if __name__ == '__main__':

    #the index pairs represent the starting index and end index of the rows to parse in the csv for each process.

    #index_pairs = [[0,922], [922,1844], [1844,2766], [2766,3689]]
    #index_pairs = [[200,250],[1000,1050], [2350,2400], [3050,3100]]
    index_pairs = [100,300]

    # multiprocessing option (index pairs must be 2D list)
    #pool = mp.Pool(mp.cpu_count())
    #pool.map(statmorphWrapper, index_pairs)
    
    # single process option (index pairs must be 1D list)
    statmorphWrapper(index_pairs)

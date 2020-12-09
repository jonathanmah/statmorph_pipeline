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
    size = 30*size
    x = (i-xCenter)
    y = (j-yCenter)
    a = math.sqrt((size**2)/ratio)
    b = (size**2)/a

    if(  (((((x*math.cos(math.radians(angle))) - (y*math.sin(math.radians(angle))))**2)  / (a**2)) + ((((y*math.cos(math.radians(angle))) + (x*math.sin(math.radians(angle))))**2)  / (b**2))) <= 1):
        return True


def statmorphWrapper(index_pairs):

    #sys.stdout = open('/home/ubuntu/statmorph_pipeline/output.txt', 'a')

    # list of too big galaxies (DELETE AFTER) temp remove  NGVSJ12:29:46.77+08:00:01.8
    toobig = ['NGVSJ12:28:03.88+09:48:13.4', 'NGVSJ12:28:43.31+11:45:18.1', 'NGVSJ12:28:57.56+13:14:31.0', 'NGVSJ12:29:00.04+13:58:42.5', 'NGVSJ12:29:21.29+08:09:23.8',
        'NGVSJ12:29:48.87+13:25:46.0', 'NGVSJ12:29:53.55+14:04:06.9', 'NGVSJ12:30:10.33+10:46:46.2', 'NGVSJ12:30:49.42+12:23:28.0',
        'NGVSJ12:32:14.22+10:15:05.3', 'NGVSJ12:34:02.98+07:41:57.5', 'NGVSJ12:34:06.08+11:19:16.5', 'NGVSJ12:35:30.58+12:13:14.7', 'NGVSJ12:35:37.95+12:15:50.4',
        'NGVSJ12:35:39.81+12:33:22.9', 'NGVSJ12:36:26.98+11:26:21.2', 'NGVSJ12:36:53.38+07:14:47.7', 'NGVSJ12:36:54.86+12:31:12.5', 'NGVSJ12:37:30.56+09:33:18.4',
        'NGVSJ12:41:32.75+07:18:53.6', 'NGVSJ12:42:02.26+11:38:49.0', 'NGVSJ12:42:47.43+11:26:33.0', 'NGVSJ12:43:39.97+11:33:09.7', 'NGVSJ12:44:31.98+11:11:25.8',
        'NGVSJ12:52:17.50+11:18:49.9', 'NGVSJ12:52:55.97+11:13:51.5'
    ]

    toobigfast = ['NGVSJ12:34:06.08+11:19:16.5', 'NGVSJ12:52:17.50+11:18:49.9', 'NGVSJ12:28:43.31+11:45:18.1', 'NGVSJ12:29:53.55+14:04:06.9', 'NGVSJ12:35:30.58+12:13:14.7', 
    'NGVSJ12:36:54.86+12:31:12.5', 'NGVSJ12:34:02.98+07:41:57.5', 'NGVSJ12:37:30.56+09:33:18.4']

    pickleError = ['NGVSJ12:42:02.26+11:38:49.0','NGVSJ12:29:48.87+13:25:46.0']
    #
    smallfile = ['NGVSJ12:29:21.29+08:09:23.8', 'NGVSJ12:41:32.75+07:18:53.6']
    
    df = pd.read_csv('NGVSgalaxies.csv')

    #iterates through rows in csv file containing galaxy names, sky value, and nuclei centers
    for row in df.iloc[index_pairs[0]:index_pairs[1]].itertuples(index=True, name='Pandas'):
        galaxy = row.Official_name
        base = 'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/files/vault/ngvs/data/NGVS/galaxies/'
        galaxyPath = f'{galaxy}/{galaxy}_G'
        

    # CHECK IF TOO BIG (DELETE AFTER)

        if(galaxy not in toobigfast):
            continue
      #  if(galaxy not in smallfile):
       #     continue
    #

        #remove these 3 lines after
        print(f'checking {galaxy}')
        writeFile = open(f'/mnt/scratch/check/{galaxy}.txt', 'w') 
        writeFile.write('hello')
        writeFile.close()


        startcopy = time.time()
        # try to copy galaxy files
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxy}/{galaxy}_G.fits /mnt/scratch/temp_galaxy_storage')
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxyPath}_iso_model.fits /mnt/scratch/temp_galaxy_storage')
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxyPath}_galfit_model.fits /mnt/scratch/temp_galaxy_storage')
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxyPath}_mask.fits /mnt/scratch/temp_galaxy_storage')
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxyPath}_psf.fits /mnt/scratch/temp_galaxy_storage')
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxyPath}_sig.fits /mnt/scratch/temp_galaxy_storage')
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxyPath}_iso_residual.fits /mnt/scratch/temp_galaxy_storage')
        os.system(f'vcp vos:ngvs/data/NGVS/galaxies/{galaxyPath}_galfit_residual.fits /mnt/scratch/temp_galaxy_storage')
        endcopy = time.time() - startcopy

        # if one of the required files is missing, move on to next galaxy                                                                     
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

        #if(os.path.getsize(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G.fits') > 300000000):
        #    writeFile = open(f'/mnt/scratch/toobig/{galaxy}.txt', 'w')     
        #    writeFile.write(f'{galaxy}\n')
        #    writeFile.close()
        #   clearTempFiles(galaxy)
        #    continue


        if(checkCorrupt(galaxy)):
            writeFile = open(f'/mnt/scratch/corrupt/{galaxy}.txt', 'w')     
            writeFile.write(f'{galaxy}\n')
            writeFile.close()
            clearTempFiles(galaxy)
            continue

        startseg = time.time()

        hdu = fits.open(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G.fits')
        im_header = hdu[0].header
        im_data = hdu[0].data

         # sky subtraction from original image
        sky_data = np.zeros(np.shape(im_data))
        sky_data += row.SKY
        im_sky_subtracted = im_data - sky_data

        # calculate nucleus center
        ra = row.NGVS_ra
        dec = row.NGVS_dec
        mywcs = wcs.WCS(im_header)
        xCenter, yCenter = mywcs.all_world2pix([[ra, dec]], 0)[0]

        mask_data = fits.getdata(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_mask.fits')
        if path.isfile(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_model.fits'):
            original_model_data = fits.getdata(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_iso_model.fits')
        else:
            original_model_data = fits.getdata(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_galfit_model.fits')

        #If nucleus exists then mask nucleus
        if(row.Nuc_Flag == 1):
            for i in range(len(mask_data)):
                for j in range(len(mask_data)):
                    if(((i-xCenter)**2) + ((j-yCenter)**2) <= (5**2)):
                        mask_data[i][j] = 100
                    
        
        ellipse_data = np.zeros(np.shape(original_model_data))
        # calculate median of original model values within 10 pixel radius from nucleus center
        pixelList = []
        for i in range(len(original_model_data)):
            for j in range(len(original_model_data)):
                    if(((i-xCenter)**2) + ((j-yCenter)**2) <= (10**2) and mask_data[i][j] != 100):
                        pixelList.append(original_model_data[i][j])   

        median = statistics.median(pixelList)

        # Create Segmentation Map
        seg_data = np.zeros(np.shape(original_model_data))

        # If median is greater than 2*sky value then create segmentation map from original model values greater than 1.4*sky value within ellipse area
        isEmpty = True

        if(median > 2*row.SKY):
            for i in range(len(original_model_data)):   
                for j in range(len(original_model_data)):
                        if(inEllipse(i,j,xCenter,yCenter,row.Size,row.AxisRatio,row.PA) and original_model_data[i][j] > (1.4*row.SKY)):
                            seg_data[i][j] = 100
                            isEmpty = False
        # If median is less than 2*sky value then create segmentation map from original model values greater than 1.1*sky value within ellipse area
        else:
            for i in range(len(original_model_data)):
                for j in range(len(original_model_data)):
                        if(inEllipse(i,j,xCenter,yCenter,row.Size,row.AxisRatio,row.PA) and original_model_data[i][j] > (1.1*row.SKY)):
                            seg_data[i][j] = 100
                            isEmpty = False

        psf = fits.getdata(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_psf.fits')
        weightmap = fits.getdata(f'/mnt/scratch/temp_galaxy_storage/{galaxy}_G_sig.fits')
        mask_data = np.array(mask_data, dtype=bool)

        endseg = time.time() - startseg



    ###########################################

    #    plt.imsave(f'/mnt/scratch/output/{galaxy}.png',seg_data, origin="lower", cmap='Greys')
    #    os.system(f'vcp /mnt/scratch/output/{galaxy}.png vos:ngvs/data/STATMORPH/4*size_segmaps/{galaxy}.png')
    #    if path.isfile(f'/mnt/scratch/output/{galaxy}.png'):
    #        os.system(f'rm /mnt/scratch/output/{galaxy}.png')

    #    hdu.close()
    #    clearTempFiles(galaxy)

        if(isEmpty):
            writeFile = open(f'/mnt/scratch/emptyseg/{galaxy}.txt', 'w')     
            writeFile.write('empty')
            writeFile.close()
            clearTempFiles(galaxy)
            continue

     #   continue

    ############################################


        start_time = time.time()
        source_morphs = statmorph.source_morphology(im_sky_subtracted, seg_data, mask=mask_data, weightmap=weightmap, psf=psf)
        end_time = time.time() - start_time

        morph = source_morphs[0]

        startmodelcreate = time.time()

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
        ax.imshow(original_res_data, norm=normOriginalRes, cmap='gray', origin='lower')
        ax.set_title('Original Residual', fontsize=15)

        ax= plt.subplot(gs[1,2])
        ax.imshow(output_res_data, norm=normOutputRes, cmap='gray', origin='lower')
        ax.set_title('Output Residual', fontsize=15)

        endfig = time.time() - startfig


        # save galaxyname - time : morphruntime - Flag=1 - SersicFlag=0
        #fig.savefig(f'/mnt/scratch/output/{galaxy}_time:{round(morph.runtime, 2)}_Flag={morph.flag}_SersicFlag={morph.flag_sersic}.png', facecolor='w', edgecolor='w', transparent=False, bbox_inches='tight')
        
        #DECEMBER1 NAME
        #fig.savefig(f'/mnt/scratch/output/{galaxy}_time:{round(end_time, 2)}_Flag={morph.flag}_SersicFlag={morph.flag_sersic}.png', facecolor='w', edgecolor='w', transparent=False, bbox_inches='tight')
        #fig.savefig(f'/mnt/scratch/output/{galaxy}_time:{round(end_time, 2)}_size={row.Size}_median={median}_2*sky={2*row.SKY}_sky={row.SKY}.png', facecolor='w', edgecolor='w', transparent=False, bbox_inches='tight')
        fig.savefig(f'/mnt/scratch/output/{galaxy}_sourcemorph:{round(end_time, 2)}_seg={round(endseg, 2)}_REfactor=30_mag={round(row.principleg_mag_cg, 3)}.png', facecolor='w', edgecolor='w', transparent=False, bbox_inches='tight')

        plt.close(fig)
        # save as multi extension fits file
        #primary_hdu = fits.PrimaryHDU(im_sky_subtracted, header=im_header)
        #image_hdu = fits.ImageHDU(model_data)
        #image_hdu2 = fits.ImageHDU(output_res_data)

        #hdul = fits.HDUList([primary_hdu, image_hdu, image_hdu2])

        # copy file to vospace
        #hdul.writeto(f'/mnt/scratch/output/{galaxy}_output.fits', overwrite=True)
        #os.system(f'vcp /mnt/scratch/output/{galaxy}_output.fits vos:ngvs/data/STATMORPH/FITS_output/{galaxy}_output.fits')
        #if path.isfile(f'/mnt/scratch/output/{galaxy}_output.fits'):
        #    os.system(f'rm /mnt/scratch/output/{galaxy}_output.fits')

        
        #SAVE FILE AS FLAGS
#        os.system(f'vcp /mnt/scratch/output/{galaxy}_time:{round(end_time, 2)}_Flag={morph.flag}_SersicFlag={morph.flag_sersic}.png \
 #       vos:ngvs/data/STATMORPH/filesize_bug/{galaxy}_time:{round(end_time, 2)}_Flag={morph.flag}_SersicFlag={morph.flag_sersic}.png')
  #      if path.isfile(f'/mnt/scratch/output/{galaxy}_time:{round(end_time, 2)}_Flag={morph.flag}_SersicFlag={morph.flag_sersic}.png'):
   #         os.system(f'rm /mnt/scratch/output/{galaxy}_time:{round(end_time, 2)}_Flag={morph.flag}_SersicFlag={morph.flag_sersic}.png')


        #SAVE FILE AS MEDIAN/SKY VALUES
      #  os.system(f'vcp /mnt/scratch/output/{galaxy}_time:{round(end_time, 2)}_size={row.Size}_median={median}_2*sky={2*row.SKY}sky={row.SKY}.png \
       # vos:ngvs/data/STATMORPH/memory_bug/{galaxy}_time:{round(end_time, 2)}_size={row.Size}_median={median}_2*sky={2*row.SKY}sky={row.SKY}.png')
       # if path.isfile(f'/mnt/scratch/output/{galaxy}_time:{round(end_time, 2)}_size={row.Size}_median={median}_2*sky={2*row.SKY}sky={row.SKY}.png'):
       #     os.system(f'rm /mnt/scratch/output/{galaxy}_time:{round(end_time, 2)}_size={row.Size}_median={median}_2*sky={2*row.SKY}sky={row.SKY}.png')

        #SAVE FILE AS RUNNING TIMES
        os.system(f'vcp /mnt/scratch/output/{galaxy}_sourcemorph:{round(end_time, 2)}_seg={round(endseg, 2)}_REfactor=30_mag={round(row.principleg_mag_cg, 3)}.png \
        vos:ngvs/data/STATMORPH/RE_factor_output/{galaxy}_sourcemorph:{round(end_time, 2)}_seg={round(endseg, 2)}_REfactor=30_mag={round(row.principleg_mag_cg, 3)}.png')
        if path.isfile(f'/mnt/scratch/output/{galaxy}_sourcemorph:{round(end_time, 2)}_seg={round(endseg, 2)}_REfactor=30_mag={round(row.principleg_mag_cg, 3)}.png'):
            os.system(f'rm /mnt/scratch/output/{galaxy}_sourcemorph:{round(end_time, 2)}_seg={round(endseg, 2)}_REfactor=30_mag={round(row.principleg_mag_cg, 3)}.png')

        hdu.close()
        clearTempFiles(galaxy)
        print(f'complete {galaxy}')
        
if __name__ == '__main__':#27 last run
    #start_time = time.time()

    #the index pairs represent the starting index and end index of the rows to parse in the csv for each process.
    #index_pairs = [[0,922], [922,1844], [1844,2766], [2766,3689]]
    #index_pairs = [[200,300],[1000,1100], [2300,2400], [3000,3100]]
    index_pairs = [[100,250],[250,375], [375,500]]
    #index_pairs = [100,500]
    #pool = mp.Pool(mp.cpu_count())
    pool = mp.Pool(3)
    pool.map(statmorphWrapper, index_pairs)
    #pool.map(intWrapper, index_pairs)
    #statmorphWrapper(index_pairs)

    #print("--- %s seconds ---" % (time.time() - start_time))
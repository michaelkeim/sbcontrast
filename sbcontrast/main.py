# -----------------------------------------------------------------------------
# Overview
# -----------------------------------------------------------------------------

# Authors: Keim, M. A., van Dokkum, P., Li, J.

# Email: michael [dot] keim [at] yale [dot] edu

# Description: The sbcontrast code is a method to obtain surface brightness
# limits that accurately reflect an images depth at a given spatial scale.
# Rather than relying on the naive Poisson expectation, this code will 
# estimate limitations presented by large scale variations. A full description
# of this method is given in Keim et al. 2022.

# Usage: Note that genuine sources should be masked - otherwise the reported
# limit will not be as deep as the true value. Thus, proper masks must be 
# generated prior to sbcontrast. Note also that the scale may be rounded so
# that the binning factors are integers. Moreover, the accuracy of
# this limit is limited by data reduction - for scales exceeding that of 
# background subtraction, genuine features above the calculated limit may 
# have been removed. Finally, at the single pixel scale sbcontrast may 
# diverge from the true rms due to correlations introduced by re-sampling.

# -----------------------------------------------------------------------------
# Imports & Dependencies
# -----------------------------------------------------------------------------

import argparse
import numpy             as     np
from   warnings          import filterwarnings
from   astropy.stats     import biweight_midvariance
from   astropy.stats     import biweight_location
from   astropy.io        import fits
from   astropy.wcs       import WCS
from   astropy.wcs.utils import proj_plane_pixel_scales

# -----------------------------------------------------------------------------
# Function to Calculate Limit
# -----------------------------------------------------------------------------

def sblimit( image, mask, pix_scale, zeropoint, sigma=1.0, scale_arcsec=60,
				minfrac=0.8, minback=6, minrun=True, verbose=True ):
	'''
	Summary
	-------
	Method to find the surface brightness limit on a given angular scale.

	Parameters
	----------
	image (numpy 2D array): input image.
	mask (numpy 2D array): mask map; nonzero pixels are masked.
	pix_scale (float): pixel scale of the image in arcsec/pixel.
	zeropoint (float): photometric zeropoint of the input image.
	sigma (float): detection threshold for the limit.
		An N=1sigma limit is calculated by default.
	scale_arcsec (float): angular scale of the limit in arcsec.
		Should be greater than pix_scale. Limit is calculated across
		60 arcsec x 60 arcsec squares by default.
	minfrac (float): acceptable fraction of masked pixels. Bins with
		less than or equal to minfrac fraction of unmasked pixels are
		discarded. Should be less than 1, set to 80% by default.
	minback (int): acceptable fraction of undiscarded bins for
		background calculation. Bins with less than or equal to minback
		surrounding undiscarded bins are not used to find the limit.
		Should be less than 8, set to 6 by default.
	minrun (bool): whether to minimize speed (True) or RAM usage (False).
		Minimizes runtime by default.
	verbose (bool): whether print out results, set to True by default.
	'''

	# Image Binning
	scale_pix = scale_arcsec / pix_scale
	scale_x = np.array([scale_pix,    int(scale_pix), 
					int(scale_pix),   int(scale_pix)+1])
	scale_y = np.array([scale_pix,    int(scale_pix),
					int(scale_pix)+1, int(scale_pix)+1])
	area = scale_x*scale_y
	area -= area[0]
	area = abs(area)[1:] # d1, d2, d3
	bin_x = int(scale_x[np.argmin(area)+1])
	bin_y = int(scale_y[np.argmin(area)+1])
	area_ratio = bin_x*bin_y / scale_pix**2
	if verbose and np.round(area_ratio, 4) != 1.:
		print('Warning: scale must be rounded for integer binning factors.')
		print('Used bin area / True bin area = {:.4f}'.format(area_ratio))
	nbins_x = int(image.shape[1] / bin_x)
	nbins_y = int(image.shape[0] / bin_y)

	# Initialize Maps
	im_loc   = np.zeros((nbins_y, nbins_x)) # Binned Map
	im_frac  = np.zeros((nbins_y, nbins_x)) # Masked Fraction
	im_fluct = np.zeros((nbins_y, nbins_x)) # Contrast/Fluctuation Map

	# Calculate Bin Values
	mask[np.nonzero(mask)] = 1.
	for i in range(nbins_x - 1):
		for j in range(nbins_y - 1):
			x1, x2, y1, y2 = i*bin_x, (i+1)*bin_x, j*bin_y, (j+1)*bin_y
			im_sec = image[y1:y2, x1:x2]
			im_mask_sec = mask[y1:y2, x1:x2]
			im_sec_in = im_sec[(im_mask_sec == 0)]
			if im_sec_in.size > 0:
				im_loc[j, i] = biweight_location(im_sec_in)
			im_frac[j, i] = 1 - float(im_sec_in.size) / float(im_sec.size)

	# Remove bins with too many masks
	im_loc[im_frac>1.-minfrac] = np.nan


	# Find local background using supermatrix
	if minrun:

		# Manually filter overly masked background
		filterwarnings(action='ignore', 
			message='All-NaN slice encountered')
	
		# Find local background from surrounding bins
		r, c = np.indices((nbins_y-2, nbins_x-2))
		bck_loc = biweight_location([im_loc[r+(int(i/3)+2)%3, c+(i+2)%3] for i in range(8)], axis=0, ignore_nan=True)
	
		# Discard poorly determined backgrounds
		bck_loc[np.sum(~np.isnan([im_loc[r+(int(i/3)+2)%3, c+(i+2)%3] for i in range(8)]), 0) <= minback] = np.nan

		# Calculate fluctuations
		im_fluct[1:-1,1:-1] = im_loc[1:-1,1:-1] - bck_loc

	# Find local background using loop
	else:

		# Calculate local background and fluctuation
		for i in range(nbins_x-2):
			for j in range(nbins_y-2):
				bck_loc = im_loc[j:j+3,i:i+3]
				bck_loc = np.delete(bck_loc.flatten(), 4)
				im_fluct[j+1,i+1] = im_loc[j+1,i+1] - biweight_location(bck_loc, ignore_nan=True) \
					if len(bck_loc[~np.isnan(bck_loc)]) > minback else np.nan

	# Remove edge
	im_fluct[:,[0,-1]] = im_fluct[[0,-1]] = np.nan

	# Find limit
	sig_adu = np.sqrt(biweight_midvariance(im_fluct, ignore_nan=True))*sigma/np.sqrt(1.125)

	# For the standard deviation of standard deviation, see:
	# https://stats.stackexchange.com/questions/631/standard-deviation-of-standard-deviation
	dsig_adu = sig_adu / np.sqrt(2.*(im_fluct.size - 1.))

	# Convert limit in ADU to magnitudes
	sb_lim = zeropoint - 2.5*np.log10(sig_adu / pix_scale**2)
	dsb_lim = 2.5*np.log10(1.+1./np.sqrt(im_fluct.size))

	if verbose:
			print('The {0:.1f}-sigma limit on {1} arcsec scales is {2:.4f} +- {3:.4f} mag/arcsec2.'.format(sigma, scale_arcsec, sb_lim, dsb_lim))

	return (sb_lim, dsb_lim), (sig_adu, dsig_adu), [im_fluct, im_loc]

# -----------------------------------------------------------------------------
# Load Files and Run
# -----------------------------------------------------------------------------

def sbc():

	# Initialize Argument Parser
	parser = argparse.ArgumentParser(description=
		'Method to find the surface brightness limit on a given angular scale.')
	parser.add_argument('image', type=str, help=
		'Required: Fits file to calculate the depth of.')
	parser.add_argument('-masks', type=str, help=
		'Optional: Map of masks. Nonzero pixels are masked.')
	parser.add_argument('-N', default='1', type=float, help=
		'Optional: Detection threshold for the limit (1sigma by default).')
	parser.add_argument('-s', default='60.', type=float,  help=
		'Optional: Angular scale of the limit in arcsec (60 by default).')
	parser.add_argument('-pix', type=float, help=
		'Optional: Pixel scale of the image in arcsec/pixel. '+
		'Read from header if not provided.')
	parser.add_argument('-zp', type=float, help=
		'Optional: Photometric zeropoint in mag/arcsec2. '+
		'Attempts to read from header if not provided.')
	parser.add_argument('-zpref', type=str, help=
		'Optional: Keyword for zeropoint in image header.')
	parser.add_argument('-binmap', type=str, help=
		'Optional: Name to save binned array (X.npy; not saved by default).')
	parser.add_argument('-fluctmap', type=str, help=
		'Optional: Name to save fluctuation array (X.npy; not saved by default).')
	parser.add_argument('-hdu', default='0', type=int, help=
		'Optional: Index of HDU object containing image data (0 by default).')
	parser.add_argument('-hdum', default='0', type=int, help=
		'Optional: Index of HDU object containing masks (0 by default).')
	args = parser.parse_args()	

	# Open File
	hduli  = fits.open(args.image)
	header = hduli[args.hdu].header if args.hdu is not None else hduli[0].header	

	# Pixel Size
	if args.pix is not None:
		pix = args.pix
	else:
		pix = proj_plane_pixel_scales(WCS(header))[0]*3600. # Assumes square pixels	

	# Zeropoint
	if args.zp is not None:
		zp = args.zp
	elif args.zpref is not None:
		zp = header[args.zpref]
	else:
		keys = list(header.keys())
		zpcheck = False
		knownrefs = ['REFZ','PHOTZ', 'Z', 'CCDZ', 'PHOTOZ']
		knownrefs.extend([ref+'P' for ref in knownrefs])
		knownrefs.extend([ref+'T' for ref in knownrefs])
		knownrefs.append('ZEROPOINT')
		for key in knownrefs:
			if key in keys:
				zp = header[key]
				zpcheck = True
				break
		if zpcheck == False:
			print(knownrefs)
			raise TypeError('Error obtaining zeropoint from header. '+
							'Please enter value manually with -zp '+
							'or provide correct reference with -zpref.')	

	# Get Data
	datai  = hduli[args.hdu].data if args.hdu is not None else hduli[0].data
	if args.masks is not None:
		hdulm  = fits.open(args.masks)
		datam  = hdulm[args.hdum].data if args.hdum is not None else hdulm[0].data
		hdulm.close()
	else:
		datam  = np.zeros(datai.shape)	

	# Close File
	hduli.close()	

	# Run
	temp = sblimit(datai, datam, pix, zp, sigma=args.N, scale_arcsec=args.s)	

	# Save Files
	if args.fluctmap is not None:
		np.save(args.fluctmap+'.npy', temp[2][0])
	if args.binmap is not None:
		np.save(args.binmap+'.npy', temp[2][1])

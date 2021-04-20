"""
==========================================
SeqFISH Image Alignment and Quantification
==========================================

Sub-pixel localization of SeqFISH data 
Image Alignment using fluorescent beads as a fiducial marker

"""

import os
import argparse
import numpy as np
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt
import imageio
#from skimage.feature import peak_local_max
from scipy import spatial
from skimage import data
from skimage.feature import register_translation
from skimage.feature.register_translation import _upsampled_dft
from scipy.ndimage import fourier_shift
from scipy.ndimage import gaussian_filter
from photutils import DAOStarFinder
import math

import time

def nearest_neighbour(points_a, points_b):
    tree = spatial.cKDTree(points_b)
    return tree.query(points_a)

def smDot(rawImg, gsigma, threshold):
	gImg = gaussian_filter(rawImg, sigma = gsigma)
	daofind = DAOStarFinder(threshold, 3)
	source = daofind(gImg)
	if (source is not None):
		if(len(source)):
			stars = np.transpose((source['xcentroid'], source['ycentroid']))
		else:
			stars = np.array([[0,0], [10000,10000]])
	else:
		stars = np.array([[0,0], [10000,10000]])
	
	return stars
	
def AlignDot(mCoordi, wCoordi, shift, thres_dis):
	wCoordi = np.add(wCoordi, shift)
	distance, index = nearest_neighbour(mCoordi, wCoordi)
	return np.count_nonzero(distance < thres_dis)

def mRegTrans(mCoordi, wCoordi, thres_dis):
	tree = spatial.cKDTree(wCoordi)
	d, i = tree.query(mCoordi, k = 100, distance_upper_bound = 100)

	dxdy = np.array([[0,0,0], [0,0,0]])
	for x in range(0,mCoordi.shape[0]):
		list = i[x]
		#print('before list:', list)
		list = list[list<wCoordi.shape[0]]
		#print('after list:', list)
		for y in range(0, len(list)):
			tdxdy = np.subtract(mCoordi[x], wCoordi[list[y]])
			shift = np.add(wCoordi, tdxdy)
			#print('shift: ', len(shift), tdxdy)
			distance, index = nearest_neighbour(mCoordi, shift)
			overlap = np.count_nonzero(distance < thres_dis)
			#print('check: ', distance, overlap)
			tdxdy = np.hstack((tdxdy, overlap))
			#print('tdxdy: ', tdxdy)
			dxdy = np.vstack((dxdy, tdxdy))

	print ('test', dxdy.shape[1])
	maxIndex = np.argmax(dxdy[:,2])
	print('Final: ', dxdy.shape[1], dxdy[maxIndex,0:2])
	return dxdy[maxIndex,0:2], dxdy[maxIndex,2]

# initiate the parser
parser = argparse.ArgumentParser()

# add long and short argument
parser.add_argument("--dirRaw", "-R", help="directory of raw images")
parser.add_argument("--dirLocal", "-L", help="directory of analysis results")
parser.add_argument("--startFOV", "-a", help="Starting Field of View")
parser.add_argument("--lastFOV", "-z", help="Last Field of View")
parser.add_argument("--m6ACyc", "-m", help="Number of m6A Cycles")
#parser.add_argument("--InoCyc", "-I", help="Number of Inosine Cycles")
parser.add_argument("--Cycles", "-c", help="Number of SeqFISH Cycles")
parser.add_argument("--thresSeq", "-s", help="threshold of seq intensity signal")
parser.add_argument("--thresDis", "-d", help="threshold of colocalization distance")
parser.add_argument("--img", "-i", help="image types [spe, tif]")
parser.add_argument("--fiducial", "-f", help="use fiducial markers")

# read arguments from the command line
args = parser.parse_args()


if args.dirRaw:
	dirRaw = args.dirRaw
else:
	dirRaw = '/Images/SeqFISH/'

if args.dirLocal:
	dirLocal = args.dirLocal
else:
	dirLocal = '/Analysis/SeqFISH/'

if args.m6ACyc:
	m6A = int(args.m6ACyc)
else:
	m6A = 0		# m6A detection cycle number

Ino = 0				# Inosine detection cycle number

if args.Cycles:
	Cyc = int(args.Cycles)
else:
	Cyc = 0		# How many of SeqFISH cycles?

if args.startFOV:
	StartFOV = int(args.startFOV)
else:
	StartFOV = 1		# Start FOV number
	
if args.lastFOV:
	NumFOV = int(args.lastFOV)
else:
	NumFOV = 5			# How many of FOVs?

TF_Ab = 10			# Noise Tolerance_Red for Ab
TF_R0 = 7			# Noise Tolerance_Red for Template
if args.thresSeq:
	TF_R1 = int(args.thresSeq)			# Noise Tolerance_Red_SeqFISH
else:
	TF_R1 = 3			# Noise Tolerance_Red_SeqFISH

TF_M = 100			# Noise Tolerance_FluoSphere Beads

if args.thresDis:
	thres_dis = int(args.thresDis)
else:
	thres_dis = 3		# Pixel distance for co-localization

if args.img:
	Img_type = args.img
else:
	Img_type = 'tif'	# Image type [spe, tif]

if args.fiducial:
	FM = int(args.fiducial)
else:
	FM = 1	# Using Fiducial Marker [TRUE, FALSE]

GF_sigma = 1


print('Load Dir: ', dirRaw)
print('Save Dir: ', dirLocal)
print('m6A cycles: ', m6A)
print('Ino cycles: ', Ino)
print('SeqFISH Cycles: ', Cyc)
print('StartFOV: ', StartFOV)
print('NumFOV: ', NumFOV)
print('TF_Ab: ', TF_Ab)
print('TF_Template: ', TF_R0)
print('TF_SeqFISH: ', TF_R1)
print('TF_Fiducial Marker: ', TF_M)
print('Image Type: ', Img_type)
print('Colocalization_Pixels: ', thres_dis)
print('Gaussian Sigma: ', GF_sigma)

start_tot = time.time()

for i in range(StartFOV, NumFOV+1):
	start = time.time()

	loadDir = dirRaw + 'Position ' + str(i) + '/';
	saveDir = dirLocal + 'Position ' + str(i) + '/';
	os.makedirs(saveDir);
	
#	templateCyc = m6A + Ino;
	templateCyc = 0;
	pathSeq = str(templateCyc) + '_' + str(i) + '_N_R_Beads.' + Img_type;
	TF = TF_R0;
	
	loadFile = loadDir +pathSeq;
	mImg = imageio.imread(loadFile);

# Find the single-molecule localizations on the template image
	mTable = smDot(mImg, GF_sigma, TF)
	saveFile = saveDir + '0_smLocal.csv';
	np.savetxt(fname=saveFile, delimiter=",", X=mTable, fmt='%.3f');
	
# Prepare the Maker image from the template image
	mImg_mask = mImg[:,:] < TF_M 	# Threshold value = TF_M
	mImg = mImg[:,:] - TF_M			
	mImg[mImg_mask] = 0				# Remain over than TF_M pixel values

	#mCoordi = smDot(mImg, GF_sigma, threshold = TF)	# Fiducial Marker Localizations
	mCoordi = mTable
	dxdy = np.array([0,0,1])
	
	print('Position Number: %d' %i)
	
	totCyc = Cyc + m6A + Ino;
	if (mCoordi.shape[0] < 2 and FM > 0):
		totCyc = 0
	
	for j in range(1, totCyc+1):
		start_img = time.time()
		if (j<=m6A):
			pathSeq = str(j) + '_' + str(i) + '_N_R_m6A.' + Img_type;
			TF = TF_Ab;
		elif(j<=(m6A+Ino)):
			pathSeq = str(j) + '_' + str(i) + '_N_R_Ino.' + Img_type;
			TF = TF_Ab;
		else:
			pathSeq = str(j) + '_' + str(i) + '_N_R_SeqFISH.' + Img_type;
			TF = TF_R1;
		
		loadFile = loadDir + pathSeq;
		wImg = imageio.imread(loadFile);

# Save the single-molecule localizations on the Ab and/or Seq images		
		wCoordi = smDot(wImg, GF_sigma, TF)
		saveFile = saveDir + str(j) + '_smLocal.csv';
		np.savetxt(fname=saveFile, delimiter=",", X=wCoordi, fmt='%.3f');
		
# Prepare the Maker image from the Ab and/or Seq image
		wImg_mask = wImg[:,:] < TF_M 	# Threshold value = TF_M
		wImg = wImg[:,:] - TF_M			
		wImg[wImg_mask] = 0				# Remain over than TF_M pixel values
		
		wCoordi_fm = smDot(wImg, GF_sigma, threshold = TF)	# Fiducial Marker Localizations

# Print the offset information		
		if (FM > 0):
			shift, error, diffphase = register_translation(mImg, wImg, 100)
			if (math.sqrt( ((shift[0])**2) + ((shift[1])**2)) > 100):
				shift, Overlap_dots = mRegTrans(mCoordi, wCoordi_fm, thres_dis)
			else:
				shift = shift[::-1]
				Overlap_dots = AlignDot(mCoordi, wCoordi, shift, thres_dis)
			
			tdxdy = np.hstack((shift, Overlap_dots/mCoordi.shape[0]))
			dxdy = np.vstack((dxdy, tdxdy))
		else:
			shift = np.array([0,0])
			tdxdy = np.hstack((shift, 1))
			dxdy = np.vstack((dxdy, tdxdy))

		tCoordi = mTable[:,0:2]
		wCoordi = np.add(wCoordi, shift)	# Shift by dxdy
		distance, index = nearest_neighbour(tCoordi, wCoordi)	# find the nearest dots
		mTable = np.column_stack((mTable, distance))	# Save the distances of the nearest dots
		
		distance, index = nearest_neighbour(wCoordi, tCoordi)	# find the nearest dots
		wCoordi_new = wCoordi[distance[:] >= thres_dis, :]		# distance filtering and find new spots
		wCoordi_nArray = np.full((wCoordi_new.shape[0], mTable.shape[1]), 10.00)	# empty array for new spots
		wCoordi_nArray[:,0:2] = wCoordi_new						# Save the coordinates of new spots
		wCoordi_nArray[:,mTable.shape[1]-1] = 0					# Save the nearest distance as 0, because they are it self. 
		mTable = np.row_stack((mTable, wCoordi_nArray))			# integrate the template and new spot list
		print('Image:', pathSeq, 'offsets: %.2f, %.2f, %.4f' %(tdxdy[0], tdxdy[1], tdxdy[2]), '		Spots: %d, %d, %.2f, %d' %(wCoordi.shape[0], wCoordi_new.shape[0], (wCoordi_new.shape[0]/wCoordi.shape[0])*100, mTable.shape[0]))
		
		end_img = time.time()
		print(j,'_cycle_time: ', end_img - start_img)
		
	if (totCyc > 0):	
	# Save the offset information		
		saveFile = saveDir + 'dxdy.csv'
		np.savetxt(fname = saveFile, delimiter = ",", X = dxdy, fmt = '%f')

	# Save Location and Intensity Information
		mTable[:, 2:totCyc+2] = mTable[:, 2:totCyc+2] < thres_dis	# distance filtering
	
		saveFile = saveDir + str(i) + '_MasterDot.csv'
		np.savetxt(fname = saveFile, delimiter = ",", X = mTable, fmt = '%.2f')
	
		tTable = mTable[:,0:2]
		if (m6A>0):
			rSum = np.sum(mTable[:,2:(2+m6A)], axis=1)				# rowSum of m6A Ab.
			rSum = rSum[:] >= 1
			tTable = np.column_stack((tTable, rSum))
			#print('head:', tTable[:10,:])
	
		if (Ino>0):
			rSum = np.sum(mTable[:,(2+m6A):(2+m6A+Ino)], axis=1)	# rowSum of Ino Ab.
			rSum = rSum[:] >= 1
			tTable = np.column_stack((tTable, rSum))
	
		tTable = np.hstack((tTable, mTable[:,(2+m6A+Ino):]))
		saveFile = dirLocal + str(i) + '_SeqDot.csv'
		np.savetxt(fname = saveFile, delimiter = ",", X = tTable, fmt = '%.2f')
	
	end = time.time()
	print(i,'_FOV_time: ', end - start)
		
end_tot = time.time()		
print('Total_Running_Time: ', end_tot - start_tot)





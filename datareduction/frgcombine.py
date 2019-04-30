from astropy.io import fits
import os
import numpy as np
'''
This script is used to subtract fringe from the object image
Written by SW Zhang
2019.04.30
'''

class Fringe(object):
	'''
	this class represent fringe image, it includes combined fringe images and its scales
	'''
	def __init__(self,masklist=[],ccdlist=[]):	
		self.masklist=masklist
		self.ccdlist=ccdlist
		self.imgdic=self._fringegenerator()
        self.scaledic=self._scalegenerator()
	def _fringegenerator(self):
		'''
		generate combined fringe images
		:return: dictionary of 2D arrays of combined fringe images
		'''
		fringedic=Combine(self.masklist)
		return fringedic 
	
	def _scalegenerator(self):
		'''
		generate scale for each combined fringe images
		:return: dictionary contain scale for each combined image
		'''
		bkgdic=Combine(self.ccdlist)
		for key in bkgdic.keys():
			bkgdic[key]=np.median(bkgdic[key])
		return bkgdic

	def fitsgenerator(self):
		for key in self.imgdic.keys():
			fringe_fits=fits.HDUList()
			fringe_fits.append(fits.ImageHDU(self.imgdic[key]))
			fringe_fits.writeto(key+'fringe.fits',overwrite=True)
			print('finish '+key+'fringe.fits')
		return None


def MaskHSource(ccdproc_list):
	'''
	mask the bright and large source in image
	:param ccdproc_list: name list of object images
	:return: name list of object masked images
	'''
	mask=[]
	for name in ccdproc_list:
		sex='sex '+name+' -c jiang001.sex -PARAMETERS_NAME jiang001.param -FILTER Y -FILTER_NAME gauss_5.0_9x9.conv -SATUR_LEVEL 60000 -DETECT_MINAREA 5 -DETECT_MAXAREA 5000 -DETECT_THRESH 2.2 -ANALYSIS_THRESH 0.2 -CHECKIMAGE_TYPE -OBJECTS -CHECKIMAGE_NAME '+'obj_'+name
		os.system(sex)
		mask.append('obj_'+name)
	return mask


def Combine(imglist):
	'''
	combine image of the same chip
	:param imglist: name list images of which you would combine
	:return: dictionary contains combined image for each chip
	'''
	imgdic={'1':[],'2':[],'3':[],'4':[]}
	for name in imglist:
		try:
			img_fits=fits.open(name)
		except IOError:
			continue
		if name[-6]=='1':
			imgdic['1'].append(img_fits[0].data)
		elif name[-6]=='2':
			imgdic['2'].append(img_fits[0].data)
		elif name[-6]=='3':
			imgdic['3'].append(img_fits[0].data)
		else:
			imgdic['4'].append(img_fits[0].data)
	for key in imgdic.keys():
		img=np.median(np.array(imgdic[key]),axis=0)
		imgdic[key]=img
	return imgdic

def FitsSelect(current_list):
	'''
	select object masked fits
	:param current_list: all file name in current dictionary
	:return: name list of object masked images
	'''
	mask_list=[]
	for name in current_list:
		if 'obj' in name:
			mask_list.append(name)
	return mask_list


def SubtractFringe(ccdproc_list,mask_list):
	'''
	subtract fringe from object images
	:param ccdproc_list: name list images of which you want to remove fringe from
	:param mask_list: name list of object masked images
	:return: None
	'''
	fringe=Fringe(mask_list,ccdproc_list)
	for name in ccdproc_list:
		print('process '+name+'.'*10)
		image_fits=fits.open(name,'update')
		i=name[-6]		
		print(image_fits[0].data[10,10])
		fringemedian=fringe.scale[i]
		sciimgmedian=np.median(image_fits[0].data)
		scale=sciimgmedian/fringemedian#calculate the scale of each object images
		print(sciimgmedian,fringemedian)
		print('the scale of '+name+' is '+str(scale))
		image_fits[0].data=image_fits[0].data-scale*fringe.imgdic[i]#subtract fringe from object images
		print(image_fits[0].data[10,10])
		image_fits.close()	
	return None

def Run():
	'''
	run the whole process
	:return: None
	'''
	os.chdir('/home/zsw666/lbtpro/data/z-SLOAN/obj/ccdproc/f1542')
	currentlist=os.listdir(os.getcwd())
	ccdlist=[name for name in currentlist if '.fits' in name and 'LBT' in name and 'obj' not in name]
	#masklist=FitsSelect(MaskHSource(ccdlist))
	masklist=[name for name in currentlist if 'obj_' in name and '.fits' in name]
	SubtractFringe(ccdlist,masklist)
	return None
Run()

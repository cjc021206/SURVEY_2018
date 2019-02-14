

import os, sys
import numpy as N
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import ICRS, FK5 
import cPickle as CP

D = os.path.abspath( sys.argv[0] )
HOME = os.path.dirname( os.path.dirname( os.path.dirname( D ) ))
IMG_HOME  = HOME + '/image/'
TMPL_HOME = HOME + '/sdss/'
OP_HOME   = HOME + '/op/'
CAT_HOME  = HOME + '/cat/'
RES_HOME  = HOME + '/result/'

IMG_ARC_DIR = IMG_HOME + 'archive/'
IMG_SUB_DIR = IMG_HOME + 'sub/'
CAT_REF_DIR = CAT_HOME + 'ref/'
CAT_IMG_DIR = CAT_HOME + 'img/'
CAT_SUB_DIR = CAT_HOME + 'sub/'
CONF_HOME = os.path.dirname( os.path.dirname( D ))+ '/config/'

N_SZ    = 0
N_BAND_IMG = 2
N_BAND_REF = 6

RADIUS_SZ  = 0.3**2
RATIO_MM   = 1.5
NSIM_CENT  = 2
NSIM_BG_0  = 4
NSIM_BG_1  = 6

NSIM_BG_TH = -3

MAG_ZERO   = 22.45

N_NONZERO    = 1
BG_ABS_FALSE = 100


XYSEP_STAR  = 5.0
XYSEP_LSTAR = 20.0
MAG_LSTAR   = 12.

NAME_REF   = 'SDSS-R9'
MARK_TMIMG = 'frame'

SUF_FITS   = '.fits'
SUF_CAT    = '.cat'
SUF_TMPL   = '_tmpl'
SUF_SUB    = '_diff'
SUF_CAND   = '_cand'
NORM       = 'i'

CAT_TMPL  = "TMPL_LIST.txt"

CDS        = '/Users/cc/lib_cc/cdsclient-3.84/'
SEX_PARA   = " -c %sbatc.sex -PARAMETERS_NAME %sbatc.param -FILTER_NAME %sbatc.conv -STARNNW_NAME %sbatc.nnw -MAG_ZEROPOINT %f" %(CONF_HOME, CONF_HOME, CONF_HOME, CONF_HOME, 22.5 )
SCAMP_PARA = " -c %sscamp.conf -CDSCLIENT_EXEC %saclient " %( CONF_HOME, CDS )
MISS_PARA  = " -c %smissfits_copy_new_wcs.conf " %( CONF_HOME )
SWARP_PARA = " -c %sswarp.conf " %( CONF_HOME )
#print SEX_PARA

def init_tmpl_coor( out_name ):
	# record the ra and dec for every sdss images into a cat file
	img_tmpl_arr = [ n_temp for n_temp in os.listdir( TMPL_HOME ) if n_temp[ :len(MARK_TMIMG) ]== MARK_TMIMG ]
	img_tmpl_arr.sort()
	out_list = ''
	out_dict = {}
	for id_tmpl in img_tmpl_arr :
		id_fits = TMPL_HOME + id_tmpl
		with fits.open( id_fits ) as f_fits:
			header_fits = f_fits[0].header
		ra_fits  = header_fits['CRVAL1']
		dec_fits = header_fits['CRVAL2']
		out_list = out_list + "%s\t%8.4f\t%7.4f\n" %(id_tmpl, ra_fits, dec_fits)
	with open( out_name, 'w' ) as out_file:
		out_file.write( out_list )
	return 0

def get_tmpl_coor( in_name ):
	# load the ra and dec information of sdss template images
	if not os.path.isfile( in_name ):
		print "NO TEMPLATE IMAGE CATALOG"
		n_temp = init_tmpl_coor( in_name )
	id_tmpl_arr, ra_tmpl_arr, dec_tmpl_arr = N.loadtxt( in_name, dtype=str, unpack=True )
	#print id_tmpl_arr
	band_tmpl_arr = N.array([ n_temp[ N_BAND_REF ] for n_temp in id_tmpl_arr ])
	ra_tmpl_arr  = ra_tmpl_arr.astype( N.float )
	dec_tmpl_arr = dec_tmpl_arr.astype( N.float )
	return [ id_tmpl_arr, band_tmpl_arr, ra_tmpl_arr, dec_tmpl_arr ]

def init_ref_cat( name_img_0, name_cat ):
	name_img    = OP_HOME + name_img_0
	name_base   = name_img_0.split('.')[0]
	name_cat_01 = OP_HOME + name_base + SUF_CAT
	band_img    = name_img_0.split('_')[ N_BAND_IMG ]
	comd_sex_01 = 'sex '+ name_img +' -CATALOG_NAME '+ name_cat_01 + SEX_PARA
	comd_scamp_01 = 'scamp '+ name_cat_01 +' -ASTREF_CATALOG SDSS-R9 -ASTREF_BAND '+ band_img + SCAMP_PARA +' -SAVE_REFCATALOG Y -REFOUT_CATPATH '+OP_HOME

	os.system( comd_sex_01 )
	os.system( comd_scamp_01 )
	SZ_sdss_cat   = [ n_temp for n_temp in os.listdir( OP_HOME ) if n_temp[ :len(NAME_REF) ]== NAME_REF ]
	mv_comd       = 'mv '+OP_HOME + SZ_sdss_cat[0] +' '+ name_cat
	os.system( mv_comd )
	return 0;

def init_wcs( name_img_0 ):
	name_img    = OP_HOME + name_img_0
	name_base   = name_img_0.split('.')[0]
	id_SZ       = name_base.split('_')[ N_SZ ]
	band_img    = name_base.split('_')[ N_BAND_IMG ]
	name_SZ_cat = CAT_REF_DIR + id_SZ +'_'+ band_img + SUF_CAT
	name_cat_01 = OP_HOME + name_base + SUF_CAT
	name_header = OP_HOME + name_base + '.head'
	name_tmpl   = OP_HOME + name_base + SUF_TMPL + '.head'
	#print name_img, name_SZ_cat, band_img; quit()
	if not os.path.isfile( name_SZ_cat ):
		print "NO LOCAL REFERENCE CATALOG for %s!" %(name_img_0)
		n_temp  = init_ref_cat( name_img_0, name_SZ_cat )
	comd_sex_01 = 'sex '+ name_img +' -CATALOG_NAME '+ name_cat_01 + SEX_PARA
	comd_scamp_01 = "scamp %s -ASTREF_CATALOG FILE %s -ASTREFCAT_NAME %s" %( name_cat_01, SCAMP_PARA, name_SZ_cat )
	#print comd_scamp_01; quit()
	comd_miss_01  = "missfits "+ name_img + MISS_PARA 
	#os.system()
	with fits.open( name_img ) as fits_img:
		head_img_0 = fits_img[0].header
	comd_sed_01    = "sed -i \'1i\\NAXIS1  = %i\\nNAXIS2  = %i\' %s" %( head_img_0['NAXIS1'], head_img_0['NAXIS2'], name_header )
	comd_mv_01     = "mv %s %s" %( name_header, name_tmpl )
	os.system( comd_sex_01 )
	os.system( comd_scamp_01 )
	os.system( comd_miss_01 )
	os.system( comd_sex_01 )
	os.system( comd_sed_01 )
	os.system( comd_mv_01 )
	#return name_SZ_cat
	return 0;

def combine_tmpl( name_img_0, info_tmpl_arr ):
	id_tmpl_arr, band_tmpl_arr, ra_tmpl_arr, dec_tmpl_arr = info_tmpl_arr
	name_img    = OP_HOME + name_img_0
	name_base   = name_img_0.split('.')[0]
	band_img    = name_base.split('_')[ N_BAND_IMG ]
	name_tmpl_0 = name_base + SUF_TMPL 
	with fits.open( name_img ) as fits_img:
		head_img_0 = fits_img[0].header
	ra_img  = head_img_0['CRVAL1']
	dec_img = head_img_0['CRVAL2']
	sep_00 = (ra_tmpl_arr - ra_img)**2 + ( dec_tmpl_arr - dec_img )**2
	keep_tmpl_m0 = ( sep_00< RADIUS_SZ )*( band_tmpl_arr==band_img )
	comd0 = 'swarp '
	for id_02 in id_tmpl_arr[ keep_tmpl_m0 ]:
		comd0 = comd0 + TMPL_HOME + id_02 + '[0] '
	comd_swarp = comd0 + SWARP_PARA + '-IMAGEOUT_NAME ' + OP_HOME + name_tmpl_0 + SUF_FITS
	comd_sex   = 'sex '+ OP_HOME + name_tmpl_0 + SUF_FITS +' -CATALOG_NAME '+ OP_HOME + name_tmpl_0 + SUF_CAT + SEX_PARA
	os.system( comd_swarp )
	os.system( comd_sex   )
	return 0;

def sub( name_img_0 ):
	name_base   = name_img_0.split('.')[0]
	name_img    = OP_HOME + name_img_0
	name_tmpl   = OP_HOME + name_base + SUF_TMPL +SUF_FITS
	name_outim  = OP_HOME + name_base + SUF_SUB  +SUF_FITS
	name_outcat = OP_HOME + name_base + SUF_SUB  + SUF_CAT
	comd_sub_01 = "/Users/cc/soft/hotpants/hotpants -c t -n %s -nrx 1 -nry 1 -iu 30000 -tl -100 -inim %s -tmplim %s -outim %s" %( NORM, name_img, name_tmpl, name_outim )
	comd_sex_01 = 'sex '+ name_outim +' -CATALOG_NAME '+ name_outcat + SEX_PARA
	os.system( comd_sub_01 )
	os.system( comd_sex_01 )

def comb_res( in_data, np_keep=5000 ):
	if N.shape( in_data )[0] == 1:
		val_arr = in_data[0]
		err_arr = N.ones( N.shape( in_data )[1] )
	else:
		val_arr, err_arr = in_data[:2]

	if len( val_arr )>np_keep:
		val_0 = N.median( val_arr )
		val_d = N.abs( val_arr - val_0 )
		val_order = N.argsort( val_d )
		val_arr   = val_arr[val_order][:np_keep]
		err_arr   = err_arr[val_order][:np_keep]
	power_arr = 1./ err_arr**2 
	power_sum = N.sum( power_arr )
	avg_val = N.sum( val_arr * power_arr ) / power_sum
	avg_err = ( N.sum( ( val_arr - avg_val )**2 * power_arr ) / power_sum )**0.5
	return [ avg_val, avg_err ]

def coor_reg( mag_list0, mag_list1, resd=0.0001 ):
	ra_t, dec_t, mag_t, magerr_t = mag_list0
	ra_i, dec_i, mag_i, magerr_i = mag_list1
	#print ra_t, ra_i; quit()
	order_t = N.argsort( ra_t )
	order_i = N.argsort( ra_i )
	#quit()
	f_mag  = 100.
	ra_k   = ra_t.copy()
	dec_k  = dec_t.copy()
	mag_k  = mag_t.copy()
	mags_k = N.ones( len(ra_t) ) * f_mag
	magse_k= N.ones( len(ra_t) ) * f_mag

	match_id  = -N.ones( len(ra_i) )
	match_sep = N.ones( len(ra_i) )

	j0_beg  = 0
	for i in range( len(ra_t) ):
		ra_t0   =     ra_t[ order_t[i] ]
		dec_t0  =    dec_t[ order_t[i] ]
		mag_t0  =    mag_t[ order_t[i] ]
		mage_t0 = magerr_t[ order_t[i] ]
		j1_beg  = 0

		for j in range( j0_beg, len(ra_i) ):
			sep_lit = match_sep[ order_i[j] ]
			ra_t1   = ra_i[ order_i[j] ]
			if ra_t1 > ( ra_t0 + resd ): break
			elif ra_t1 < ( ra_t0 - resd ): j0_beg = j
			dec_t1  = dec_i[ order_i[j] ]
			mag_t1  = mag_i[ order_i[j] ]
			#if abs( mag_t1 - mag_t0 )>50: continue
			mage_t1 = magerr_i[ order_i[j] ]
			if N.abs( ra_t1 - ra_t0 )>resd or N.abs( dec_t1 - dec_t0 )>resd: continue
			sep_t  = (ra_t1 - ra_t0)**2 + ( dec_t1 - dec_t0 )**2
			#print 'separation', sep_t, N.abs( ra_t1 - ra_t0 )>resd, N.abs( dec_t1 - dec_t0 )>resd
			if sep_t < sep_lit:
				mags_k[i] = mag_t0 - mag_t1
				magse_k[i]= (mage_t0**2 + mage_t1**2)**0.5 
				match_sep[ order_i[j] ] = sep_t
				match_id[ order_i[j] ]  = order_t[i]
	m_keep = mags_k < f_mag
	mag_sh, magerr_sh = comb_res( [ mags_k[m_keep], magse_k[m_keep] ], NUM_SHIFT )
	return [ mag_sh, magerr_sh, match_id ]

def SN_identif( name_img_0 ):
	name_base   = name_img_0.split('.')[0]
	name_img    = OP_HOME + name_img_0
	name_cat_inim = OP_HOME + name_base + SUF_CAT
	name_cat_tmpl = OP_HOME + name_base + SUF_TMPL + SUF_CAT
	name_cat_sub  = OP_HOME + name_base + SUF_SUB  + SUF_CAT
	name_img_sub  = OP_HOME + name_base + SUF_SUB  + SUF_FITS
	name_res_cat  = OP_HOME + name_base + SUF_CAND + SUF_CAT

	id_SZ       = name_base.split('_')[ N_SZ ]
	band_img    = name_base.split('_')[ N_BAND_IMG ]
	name_SZ_cat = CAT_REF_DIR + id_SZ +'_'+ band_img + SUF_CAT

	# read inim sex catalog
	table_inim = Table.read( name_cat_inim, hdu=2 )
	ra_inim    = N.array( table_inim['ALPHA_J2000'] )
	dec_inim   = N.array( table_inim['DELTA_J2000'] )
	mag_inim   = N.array( table_inim['MAG_AUTO'] ) 
	mage_inim  = N.array( table_inim['MAGERR_AUTO'] )
	fwhm_inim  = N.array( table_inim['FWHM_IMAGE'] ) 
	fr_inim    = N.array( table_inim['FLUX_RADIUS'] ) 
	e_inim     = N.array( table_inim['ELLIPTICITY'] ) 
	mu_inim    = N.array( table_inim['MU_MAX'] ) 
	
	M_fwhm_inim = N.median( fwhm_inim )
	E_fwhm_inim = N.median( N.abs( fwhm_inim  - M_fwhm_inim ) )
	
	M_fr_inim = N.median( fr_inim )
	E_fr = N.median( N.abs( fr_inim - M_fr_inim ) )
	
	idx_fwhm_inim = ( N.abs( fwhm_inim - M_fwhm_inim)<2.5*E_fwhm_inim ) * ( N.abs(fr_inim-M_fr_inim)<2.5*E_fr )
	
	M_e_inim = N.median(e_inim[idx_fwhm_inim])
	E_e_inim = N.median( N.abs( e_inim[idx_fwhm_inim] - M_e_inim ) )
	
	# read tmpl sex catalog
	table_tmpl = Table.read( name_cat_tmpl, hdu=2 )
	ra_tmpl    = N.array( table_tmpl['ALPHA_J2000'] )
	dec_tmpl   = N.array( table_tmpl['DELTA_J2000'] )
	mag_tmpl   = N.array( table_tmpl['MAG_AUTO'] )
	mage_tmpl  = N.array( table_tmpl['MAGERR_AUTO'] )
	x_tmpl     = N.array( table_tmpl['X_IMAGE'] )
	y_tmpl     = N.array( table_tmpl['Y_IMAGE'] )
	e_tmpl     = N.array( table_tmpl['ELLIPTICITY'] )
	
	idx_e_tmpl = ( e_tmpl<0.2 )
	
	x_tmpl   = x_tmpl[idx_e_tmpl].copy()
	y_tmpl   = y_tmpl[idx_e_tmpl].copy()
	mag_tmpl = mag_tmpl[idx_e_tmpl].copy()
	
	# read subtracted sex catalog
	table_sub = Table.read( name_cat_sub, hdu=2 )
	mag_sub   = N.array( table_sub['MAG_AUTO'] )
	mage_sub  = N.array( table_sub['MAGERR_AUTO'] )
	fwhm_sub  = N.array( table_sub['FWHM_IMAGE'] )
	e_sub     = N.array( table_sub['ELLIPTICITY'] )
	x_sub     = N.array( table_sub['X_IMAGE'] )
	y_sub     = N.array( table_sub['Y_IMAGE'] )
	ra_sub    = N.array( table_sub['ALPHA_J2000'] )
	dec_sub   = N.array( table_sub['DELTA_J2000'] )
	mu_sub    = N.array( table_sub['MU_MAX'] )
	fr_sub    = N.array( table_sub['FLUX_RADIUS'] )
	th_sub    = N.array( table_sub['THRESHOLD'] )
	
	idx_e_sub = ( e_sub<min( (M_e_inim+4.0*E_e_inim), 0.3 ) )
	l1 = len( N.nonzero(idx_e_sub)[0] )
	
	idx_fwhm_sub = ( (fwhm_sub < (M_fwhm_inim+3.5*E_fwhm_inim ) ) & ( fwhm_sub >(M_fwhm_inim-2.5*E_fwhm_inim) ) )
	l2 = len( N.nonzero(idx_fwhm_sub)[0] )
	
	idx_fr_sub = ( (fr_sub< (M_fr_inim+3.5*E_fr ) ) & ( fr_sub>(M_fr_inim-2.5*E_fr) ) )
	#idx_fr_sub = ( N.abs(fr_sub-M_fr_inim)<3.5*E_fr )
	l3 = len( N.nonzero(idx_fr_sub)[0] )
	
	idx_th_sub = th_sub>0.
	'''
	l12 = len( N.nonzero(idx_e_sub*idx_fwhm_sub)[0] )
	l13 = len( N.nonzero(idx_e_sub*idx_fr_sub)[0] )
	l23 = len( N.nonzero(idx_fr_sub*idx_fwhm_sub)[0] )
	'''
	idx_sub = idx_e_sub * idx_fwhm_sub * idx_fr_sub * idx_th_sub
	
	O_x = x_sub[idx_sub].copy()
	O_y = y_sub[idx_sub].copy()
	O_ra = ra_sub[idx_sub].copy()
	O_dec = dec_sub[idx_sub].copy()
	O_mag = mag_sub[idx_sub].copy()
	O_mage = mage_sub[idx_sub].copy()
	O_fwhm = fwhm_sub[idx_sub].copy()
	O_e = e_sub[idx_sub].copy()
	
	L = len( O_x )
	print len(x_sub), L
	flag_out = N.zeros(O_x.shape) +3

	if len(O_x)<120:
		flag_out = flag_out - 3 
		with fits.open( name_img ) as f_image:
			imdata_inim = f_image[0].data
		Y_MAX_inim, X_MAX_inim = imdata_inim.shape
		xy_sreg  = NSIM_CENT * M_fwhm_inim 
	
		xy_bgreg_0 = NSIM_BG_0 * M_fwhm_inim
		xy_bgreg_1 = NSIM_BG_1 * M_fwhm_inim
	
		id1d_x_inim = N.arange( X_MAX_inim )
		id1d_y_inim = N.arange( Y_MAX_inim )
		id2d_x_inim, id2d_y_inim = N.meshgrid( id1d_x_inim, id1d_y_inim )
		for i in range(len(O_x)):
			x_cand_0  = int(O_x[i])
			y_cand_0  = int(O_y[i])
			idx_scent_sub = ( id2d_x_inim>(x_cand_0-xy_sreg) )*( id2d_x_inim<(x_cand_0+xy_sreg) )*( id2d_y_inim>(y_cand_0-xy_sreg) )*( id2d_y_inim<(y_cand_0+xy_sreg) )
			imdata_reg    = imdata_inim[ idx_scent_sub ].copy()
			
			idx_sbg0_sub = ( id2d_x_inim>(x_cand_0- xy_bgreg_0 ) )*( id2d_x_inim<(x_cand_0+ xy_bgreg_0 ) )*( id2d_y_inim>(y_cand_0- xy_bgreg_0 ) )*( id2d_y_inim<(y_cand_0+ xy_bgreg_0 ) )
			idx_sbg1_sub = ( id2d_x_inim>(x_cand_0- xy_bgreg_1 ) )*( id2d_x_inim<(x_cand_0+ xy_bgreg_1 ) )*( id2d_y_inim>(y_cand_0- xy_bgreg_1 ) )*( id2d_y_inim<(y_cand_0+ xy_bgreg_1 ) )
			imdata_bg    = imdata_inim[ idx_sbg1_sub - idx_sbg0_sub ].copy()
			
			bg_cand_00   = N.average( imdata_bg )
			bgth_cand_00 = NSIM_BG_TH * ( N.std(imdata_bg) );
			
			zero_n = len( N.nonzero( ( imdata_reg - bg_cand_00 )< bgth_cand_00 )[0] )

			if ( zero_n> N_NONZERO or abs(bg_cand_00)> BG_ABS_FALSE ):
			#if zero_n>1:
				#print abs(t_bg), abs(av)
				flag_out[i] = 4

			idx_tstar_sub  = (N.abs( x_tmpl-O_x[i] ) < XYSEP_STAR) * (N.abs( y_tmpl-O_y[i] ) < XYSEP_STAR)
			mag_bstar_tmpl = mag_tmpl[(N.abs( x_tmpl-O_x[i] ) < XYSEP_LSTAR) * (N.abs( y_tmpl-O_y[i] ) < XYSEP_LSTAR)]
			b_i = ( mag_bstar_tmpl < MAG_LSTAR )
			if ( len( N.nonzero(idx_tstar_sub)[0]) + len( N.nonzero(b_i)[0] ) )>0 :
				flag_out[i] = 1

	idx_keep  = flag_out<4
	O_x = O_x[idx_keep].copy()
	O_y = O_y[idx_keep].copy()
	O_ra  = O_ra[idx_keep].copy()
	O_dec = O_dec[idx_keep].copy()
	O_mag  = O_mag[idx_keep].copy()
	O_mage = O_mage[idx_keep].copy()
	O_fwhm = O_fwhm[idx_keep].copy()
	O_e = O_e[idx_keep].copy()
	
	L = len( O_x )
	print '%s final sources: %i' %(name_img_0, L)

	if L>0:
		data_ref     = Table.read( name_SZ_cat, hdu=2 )
		ra_cat_ref   = data_ref['X_WORLD']
		dec_cat_ref  = data_ref['Y_WORLD']
		mag_cat_ref  = data_ref['MAG']
		merr_cat_ref = data_ref['MAGERR']
		match_data_0 = [ ra_cat_ref, dec_cat_ref, mag_cat_ref, merr_cat_ref ]
		if NORM=='i':
			match_data_1 = [ ra_inim, dec_inim, mag_inim, mage_inim ]
		if NORM=='t':
			match_data_1 = [ ra_tmpl, dec_tmpl, mag_tmpl, mage_tmpl ]
		res_match = coor_reg( match_data_0, match_data_1 )
		mag_sh_val = res_match[0]
		mag_sh_err = res_match[1]
		O_mag  = O_mag + mag_sh_val
		O_mage = ( O_mage**2 + mag_sh_err**2 )** 0.5

	out_data = [ O_x, O_y, O_ra, O_dec, O_mag, O_mage, O_fwhm, O_e, flag_out ]
	#print out_data
	#out_table = Table( out_data, names=('X_image', 'Y_image', 'RA(deg)', 'DEC(deg)', 'mag', 'Err_mag', 'fwhm(pix)', 'e', 'FLAG'), dtype=(float, float, float, float, float, float, float, float, int ) )
	#out_table .write( name_res_cat, format='ascii')
	out_line = 'X_image\tY_image\tRA(deg)\tDEC(deg)\tmag\tErr_mag\tfwhm(pix)\te\tFLAG\n'
	for j in range( L ):
		line_temp = "%5.1f\t%5.1f\t%8.4f\t%8.4f\t%6.3f\t%5.3f\t%5.2f\t%4.2f\t%i\n" %(O_x[j], O_y[j], O_ra[j], O_dec[j], O_mag[j], O_mage[j], O_fwhm[j], O_e[j], flag_out[j] )
		out_line = out_line + line_temp
	with open(name_res_cat, 'w') as res_file:
		res_file.write( out_line )
	return L

def main():
	import datetime
	info_tmpl_arr = get_tmpl_coor( CAT_HOME + CAT_TMPL )

	list_img  = [ n_temp for n_temp in os.listdir( IMG_HOME ) if n_temp[-4:]=='fits' ]
	out_line  = '\timage\tcand\ttime(total, wcs, comb, sub, det)\n'
	for name_img in list_img:
		print name_img
		t_begin = datetime.datetime.now()
		os.system( "cp %s%s %s" %( IMG_HOME, name_img, OP_HOME ) )
		init_wcs( name_img )#; quit()
		t_1 = datetime.datetime.now()
		combine_tmpl( name_img, info_tmpl_arr )
		t_2 = datetime.datetime.now()
		sub( name_img )
		t_3 = datetime.datetime.now()
		n_cand = SN_identif( name_img )
		t_finish = datetime.datetime.now()

		print name_img, 'cost time:', (t_finish - t_begin).seconds, (t_1 - t_begin).seconds, (t_2 - t_1).seconds, (t_3 - t_2).seconds, (t_finish - t_3).seconds
		new_line = "%s\t%4i\t%4i\t%4i\t%4i\t%4i\t%4i\n" %(name_img, n_cand, (t_finish - t_begin).seconds, (t_1 - t_begin).seconds, (t_2 - t_1).seconds, (t_3 - t_2).seconds, (t_finish - t_3).seconds )
		out_line = out_line + new_line
		#print t_begin, t_finish
		#quit()

	with open( 'time_record_l', 'w') as res_file:
		res_file.write( out_line )


if __name__ == '__main__':
	i = main()
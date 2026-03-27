#! /usr/bin/python
# -*- coding:utf8 -*-

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import time
import multiprocessing

fontsize=15
plt.rc("font",size=fontsize)
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
#matplotlib.use('TkAgg')

color=['#ff0000','#ff6600','#00ff00','#006600','#00ffff','#0000ff','#cc66ff','k']
fmt='os>*^hv'

clock=time.time()

epsilon=0.3
rho0=5
beta=1
LX=200
LY=200
Rd=10
rhod=25
init=0
tmax=1000
NCPU=4
multi=True
movie=False

for arg in sys.argv[1:]:
	if "-epsilon=" in arg:
		epsilon=float(arg[9:])
	elif "-rho0=" in arg:
		rho0=float(arg[6:])
	elif "-beta=" in arg:
		beta=float(arg[6:])
	elif "-LX=" in arg:
		LX=int(arg[4:])
	elif "-LY=" in arg:
		LY=int(arg[4:])
	elif "-Rd=" in arg:
		Rd=float(arg[4:])
	elif "-rhod=" in arg:
		rhod=float(arg[6:])
	elif "-init=" in arg:
		init=int(arg[6:])
	elif "-tmax=" in arg:
		tmax=int(arg[6:])
	elif "-NCPU=" in arg:
		NCPU=int(arg[6:])
	elif "-movie" in arg:
		movie=True
	else:
		print("Bad Argument: ",arg)
		sys.exit(1)
		
if NCPU==1:
	multi=False
elif NCPU<1:
	print("Bad value of NCPU: ",NCPU)
	sys.exit(1)

def timeSnap(i,j):
	if i*DT[j]<=400:
		return i*DT[j]
	else:
		return 400+(i-400/DT[j])*DT2[j]

dpi=240
EPSILON=[0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7]

if beta==1 and init==0:
	DT=[6,3,3,2,2,2,1,1,1]
	DT2=DT
	color='kkkkkkkkk'
elif beta==0.75 and init==0:
	DT=[5,3,3,2,2,2,2,2,1]
	DT2=DT
	color='kkkkkkkkk'
elif beta==1 and init==1:
	DT=[1,1,1,1,2,2,2,2,2]
	DT2=[1,1,1,1,2,2,2,3,12]
	color='wwwwwkkkk'
elif beta==2:
	DT=[3,3,3,4,4,5,7,8,15]
	DT2=DT
	color='wwwwwwwww'

figalpha='abcdefghi'

colors_map1=[(0,0,1),(1,1,1),(1,2/3.,2/3.),(1,1/3.,1/3.),(1,0,0)]
colors_map2=[(0,0,1),(1,1,1),(2/3.,1,2/3.),(1/3.,1,1/3.),(0,1,0)]

if init==0:
	cmap=matplotlib.colors.LinearSegmentedColormap.from_list('my_list1', colors_map1, N=256)
else:
	cmap=matplotlib.colors.LinearSegmentedColormap.from_list('my_list2', colors_map2, N=256)

def Snapshot(i):
	path='snapshots/figure_APM4_mag_beta=%.8g_rho0=%.8g_init=%d_%d.png'%(beta,rho0,init,i)
	if not os.path.isfile(path):
		fig=plt.figure(figsize=(9,9))
		gs=matplotlib.gridspec.GridSpec(3,3,width_ratios=[1,1,1],height_ratios=[1,1,1],left=0.06,right=0.98,bottom=0.05,top=0.97,hspace=0.25,wspace=0.25)
		
		for j,epsilon in enumerate([0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7]):
			t=timeSnap(i,j)
			ax=plt.subplot(gs[j//3,j%3])
			
			MAG=np.fromfile('data_APM4_dynamics2d/APM4_mag_beta=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_Rd=%.8g_rhod=%.8g_init=%d_t=%d.bin'%(beta,epsilon,rho0,LX,LY,Rd,rhod,init,t),dtype=np.float32).reshape(LY,LX)
			x=np.linspace(0,LX,LX)
			y=np.linspace(0,LY,LY)
			plt.pcolormesh(LX-x,y,MAG,rasterized=True,vmin=-1/3.,vmax=1,cmap=cmap)

			#plt.axis('equal')	
			plt.xlim([0,LX])
			plt.ylim([0,LY])
			plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
			plt.yticks([0,0.25*LY,0.5*LY,0.75*LY,LY])
			
			ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
			ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))

			plt.text(0.01*LX,0.04*LY,'$t=%d$'%(t),ha='left',va='center',fontsize=15,color=color[j])
			plt.text(LX,1.04*LY,'$\\beta=%.8g$, $\\epsilon=%.8g$, $\\rho_0=%.8g$'%(beta,epsilon,rho0),ha='right',va='center',fontsize=15)
			plt.text(0.99*LX,0.95*LY,'$r_d=%.8g$, $\\rho_0^d=%.8g$'%(Rd,rhod),ha='right',va='center',fontsize=15,color=color[j])
			
		plt.savefig(path,dpi=dpi)
		plt.close()		
		
		print('-snap=%d/%d -t=%.8g -tcpu=%d'%(i+1,Nsnap,t,time.time()-clock))
		del MAG,fig


os.system('mkdir -p snapshots')

i=0
t=0
ARG=[]
TIME=[]
while os.path.isfile('data_APM4_dynamics2d/APM4_mag_beta=%.8g_epsilon=%.8g_rho0=%.8g_LX=%d_LY=%d_Rd=%.8g_rhod=%.8g_init=%d_t=%d.bin'%(beta,epsilon,rho0,LX,LY,Rd,rhod,init,t)) and t<=tmax:
	ARG.append(i)
	TIME.append(t)
	i+=1
	t+=DT[0]

Nsnap=len(ARG)	
print('%d Snapshots'%Nsnap)
if multi:
	pool=multiprocessing.Pool(NCPU)
	results=pool.imap_unordered(Snapshot,ARG[::-1])
	pool.close()
	pool.join()
else:
	for i in ARG[::-1]:
		Snapshot(i)	

if movie:
	os.system('mkdir -p movies')
	os.system('ffmpeg -v quiet -stats -y -r 25/1 -i snapshots/figure_APM4_mag_beta=%.8g_rho0=%.8g_init=%d_%%01d.png -c:v h264 -r 25 -crf 25 -s %dx%d movies/movie_APM4_hydro_mag_beta=%.8g_rho0=%.8g_init=%d.mp4'%(beta,rho0,init,9*dpi,9*dpi,beta,rho0,init))

print('OK - time=%d sec'%(time.time()-clock))

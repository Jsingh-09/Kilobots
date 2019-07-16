import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import math
import os

files=[]
for (dirpath, dirnames, filenames) in os.walk('frames/'):
	files.extend(['frames/'+filenamest for filenamest in filenames if 'DS_Store' not in filenamest])
	break
files = sorted(files)
#print(files)
fig = plt.figure(figsize=(20,20))
plt.clf()

def readData(i):
	#print(i)
	plt.clf()
	data=[]
	data=np.loadtxt(files[i])
	x=zip(*data)[1]
	y=zip(*data)[2]
	COMx=np.mean(np.array(x))
	COMy=np.mean(np.array(y))
	colors=zip(*data)[0]
	size=[100*i*i for i in zip(*data)[3]]
	crossx=int(math.ceil(COMx/50.0))*50
	crossy=int(math.ceil(COMy/50.0))*50
	plt.plot([crossx,crossx],[crossy-200,crossy+200],'k--',alpha=0.5)
	plt.plot([crossx-200,crossx+200],[crossy,crossy],'k--',alpha=0.5)
	plt.plot([crossx+50,crossx+50],[crossy-200,crossy+200],'k--',alpha=0.5)
	plt.plot([crossx-200,crossx+200],[crossy+50,crossy+50],'k--',alpha=0.5)
	plt.plot([crossx-50,crossx-50],[crossy-200,crossy+200],'k--',alpha=0.5)
	plt.plot([crossx-200,crossx+200],[crossy-50,crossy-50],'k--',alpha=0.5)
	plt.plot([crossx-100,crossx-100],[crossy-200,crossy+200],'k--',alpha=0.5)
	plt.plot([crossx-200,crossx+200],[crossy-100,crossy-100],'k--',alpha=0.5)
	plt.scatter(x, y, c=colors, s=size, vmin=-0.25, vmax=2.5, edgecolors='none', cmap='spectral')
	plt.xlim([COMx-100,COMx+100])
	plt.ylim([COMy-100,COMy+100])

	plt.axis('off')
	plt.subplots_adjust(left=-0.205,bottom=-0.205,right=1.05,top=1.05)
	
anim = ani.FuncAnimation(fig,readData,frames=len(files), blit=False)
anim.save("Users/jashanbhinder/Desktop/Repulsion_first_try_nochange/video/saves/animation.mp4", fps=60, extra_args=['-vcodec', 'libx264'])


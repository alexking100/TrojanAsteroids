'''
under progress!! 
'''


from trojan import *
from matplotlib import animation
from matplotlib import rc
dxs = get_ds(asteroids)
dys = get_ds(asteroids)
solved_grids = []
xpositions = []
ypositions = []
for i, dx in enumerate(dxs):
	for j, dy in enumerate(dys):
		solved_grid = ODE(np.loadtxt(f'data/grid_{i}_{j}.txt'))
		solved_grids.append(solved_grid)
		xpositions.append(solved_grid.r_a.x)
		ypositions.append(solved_grid.r_a.y)

print(np.shape(xpositions))

# plt.plot(xpositions, ypositions)
# plt.show()

# xdata is a 2d array with shape (asteroids**2, N = 1000) 
# xdata = np.transpose(xdata)
# ydata = np.transpose(ydata)



# rc('animation', html='jshtml') # ?? what does this do

animatfig, ax = plt.subplots()


def update_animation(k):

  ax.clear()
  # xdata = xpositions[:, k]
  # ydata = ypositions[:, k]

  xdata = np.linspace(0,1,100)[k]
  ydata = np.linspace(4,5,100)[k]
  # plxdata = xpositions[k,255:257]
  # plydata = ypositions[k,255:257]
  plot=ax.plot(xdata, ydata, 'o', MarkerSize = 1, color = 'blue')
  # plot=ax.plot(plxdata, plydata, 'o', MarkerSize = 5, color = "red")

  ax.set_aspect('equal')
  # ax.set_xlim([-60,60])
  # ax.set_ylim([-60,60])
  plt.draw
  return plot
 # return (line,)


animat = animation.FuncAnimation(animatfig, update_animation, frames=np.linspace(0,1,100), interval = 20)
#animat = animation.FuncAnimation(animatfig, update_animation, frames=tvalues, interval = 20)
plt.show()
animat

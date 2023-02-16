''' 
Functions for working on AWS EC2 instance

@author: jtomf
jfarrar@whoi.edu
'''

from netrc import netrc
from subprocess import Popen
from platform import system
from getpass import getpass
import os

def setup_netrc_file():
    """PODAAC recipe for setting up netrc file for data downloads.
    Recipe obtained from PODAAC here:
    https://nasa-openscapes.github.io/earthdata-cloud-cookbook/get-started/earthdata-login.html
    """
    urs = 'urs.earthdata.nasa.gov'  # Earthdata URL endpoint for authentication
    prompts = [
        'Enter NASA Earthdata Login Username: ',
        'Enter NASA Earthdata Login Password: '
    ]

    # Determine the OS (Windows machines usually use an '_netrc' fi
    netrc_name = "_netrc" if system() == "Windows" else ".netrc"
    homeDir = os.path.expanduser("~")

    # Determine if netrc file exists, and if so, if it includes NASA Earthdata Login Credentials
    try:
        netrcDir = os.path.expanduser(f"~/{netrc_name}")
        netrc(netrcDir).authenticators(urs)[0]
        print('valid netrc file found')

    # Below, create a netrc file and prompt user for NASA Earthdata Login Username and Password
    except FileNotFoundError:
        print('netrc file not found, please login into NASA Earthdata:')
        Popen('touch {0}{2} | echo machine {1} >> {0}{2}'.format(
            homeDir + os.sep, urs, netrc_name),
              shell=True)
        Popen('echo login {} >> {}{}'.format(getpass(prompt=prompts[0]),
                                             homeDir + os.sep, netrc_name),
              shell=True)
        Popen('echo \'password {} \'>> {}{}'.format(getpass(prompt=prompts[1]),
                                                    homeDir + os.sep,
                                                    netrc_name),
              shell=True)
        # Set restrictive permissions
        Popen('chmod 0600 {0}{1}'.format(homeDir + os.sep, netrc_name),
              shell=True)

        print(f'nterc file written to {homeDir}{os.sep}{netrc_name}')

        # Determine OS and edit netrc file if it exists but is not set up for NASA Earthdata Login
    except TypeError:
        print(
            'netrc exists but is not set up for NASA Earthdata Login, please login into NASA Earthdata:'
        )
        homeDir = os.path.expanduser("~")
        Popen('echo machine {1} >> {0}{2}'.format(homeDir + os.sep, urs,
                                                  netrc_name),
              shell=True)
        Popen('echo login {} >> {}{}'.format(getpass(prompt=prompts[0]),
                                             homeDir + os.sep, netrc_name),
              shell=True)
        Popen('echo \'password {} \'>> {}{}'.format(getpass(prompt=prompts[1]),
                                                    homeDir + os.sep,
                                                    netrc_name),
              shell=True)

        print(f'nterc file written to {homeDir}{os.sep}{netrc_name}')

    return f'{homeDir}{os.sep}{netrc_name}'



def plot_ops_area(ax,**kwargs):
   """ Add polygon to show S-MODE IOP1 operations area.
         
   Inputs
   - matplotlib.pyplot.plot kwargs

   Return
   - exit code (True if OK)
   """
    # Add S-MODE IOP operations area
   '''
    New corners of polygon:
    35.790897° -125.538656°
    38.182585° -126.244555°
    38.113294° -125.438469°
    37.711709° -123.994657°
    37.795817° -123.383396°
    36.998690° -122.922778°
    '''
    
   coord = [[-125.538656,35.790897], [-126.244555,38.182585], [-125.438469,38.113294], [-123.994657,37.711709], [-123.383396,37.795817], [-122.922778, 36.998690]]
   coord.append(coord[0]) #repeat the first point to create a 'closed loop'

   xs, ys = zip(*coord) #create lists of x and y values

   if ax is None:
       ax = plt.gca()    
   # ax.plot(xs,ys,transform=ccrs.PlateCarree()) 
   ax.plot(xs,ys,**kwargs) 

   SF_lon=-(122+25/60)
   SF_lat= 37+47/60

   # mark a known place to help us geo-locate ourselves
   '''
   ax.plot(SF_lon, SF_lat, 'o', markersize=3, zorder=10, **kwargs)
   ax.text(SF_lon-5/60, SF_lat+5/60, 'San Francisco', fontsize=8, zorder=10, **kwargs)
   # ax.text(np.mean(xs)-.6, np.mean(ys)-.3, 'S-MODE ops area', fontsize=8, **kwargs)
   print(kwargs)
   '''

   return(xs,ys,ax)

def plot_ops_area_IOP2(ax,**kwargs):
   """ Add polygon to show S-MODE IOP1 operations area.
         
   Inputs
   - matplotlib.pyplot.plot kwargs

   Return
   - exit code (True if OK)
   """
    # Add S-MODE IOP1 operations area
   '''
    New corners of polygon:
    35.790897° -125.538656°
    38.182585° -126.244555°
    38.113294° -125.438469°
    37.711709° -123.994657°
    37.795817° -123.383396°
    36.998690° -122.922778°
    '''
    # Add S-MODE IOP2 operations area
   '''
    New corners of polygon:
      38.342° -126.25°
      37.707° -123.99°
      37.75° -123.354°
      37.00° -122.92°
      36.337° -124.36°
      35.00° -123.50°
      35.00° -125.353°    
    '''
    
   coord = [[38.342, -126.25], [37.707, -123.99], [37.75, -123.354], [37.00, -122.92], [36.337, -124.36], [35.00, -123.50],[35.00, -125.353]]
   coord.append(coord[0]) #repeat the first point to create a 'closed loop'

   ys, xs = zip(*coord) #create lists of x and y values

   if ax is None:
       ax = plt.gca()    
   # ax.plot(xs,ys,transform=ccrs.PlateCarree()) 
   ax.plot(xs,ys,**kwargs) 

   SF_lon=-(122+25/60)
   SF_lat= 37+47/60

   # mark a known place to help us geo-locate ourselves
   '''
   ax.plot(SF_lon, SF_lat, 'o', markersize=3, zorder=10, **kwargs)
   ax.text(SF_lon-5/60, SF_lat+5/60, 'San Francisco', fontsize=8, zorder=10, **kwargs)
   # ax.text(np.mean(xs)-.6, np.mean(ys)-.3, 'S-MODE ops area', fontsize=8, **kwargs)
   print(kwargs)
   '''

   return(xs,ys,ax)


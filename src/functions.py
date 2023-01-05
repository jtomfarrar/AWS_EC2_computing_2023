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
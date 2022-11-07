# encoding: utf-8

"""
Test module for ``ldsc.sif`` singularity build 
or ``ldsc`` dockerfile build

In case ``singularity`` is unavailable, the test function(s) should fall 
back to ``docker``.
"""

import os
import socket
import subprocess


# port used by tests
sock = socket.socket()
sock.bind(('', 0))
port = sock.getsockname()[1]

# Check that (1) singularity exist, and (2) if not, check for docker. 
# If neither are found, tests will not run.
try:
    pth = os.path.join('containers', 'ldsc.sif')
    out = subprocess.run('singularity')
    cwd = os.getcwd()
    PREFIX = f'singularity run {pth}'
    PREFIX_MOUNT = PREFIX_MOUNT = f'singularity run --home={cwd}:/home/ {pth}'
except FileNotFoundError:    
    try:
        out = subprocess.run('docker')
        pwd = os.getcwd()
        PREFIX = f'docker run -p {port}:{port} ldsc'
        PREFIX_MOUNT = (f'docker run -p {port}:{port} ' + 
            f'--mount type=bind,source={pwd},target={pwd} ldsc')
    except FileNotFoundError:
        raise FileNotFoundError('Neither `singularity` nor `docker` found in PATH. Can not run tests!')

def test_assert():
    """dummy test that should pass"""
    assert True

def test_ldsc_python():
    """test that the Python installation works"""
    call = f'{PREFIX} python --version'
    out = subprocess.run(call.split(' '))
    assert out.returncode == 0

def test_ldsc_munge():
    pwd = os.getcwd() if PREFIX.rfind('docker') >= 0 else '.'
    call = f'''{PREFIX_MOUNT} python /tools/ldsc/munge_sumstats.py --help'''
    out = subprocess.run(call.split(' '), capture_output=True)
    assert out.returncode == 0
    
def test_ldsc_ldsc():
    pwd = os.getcwd() if PREFIX.rfind('docker') >= 0 else '.'
    call = f'''{PREFIX_MOUNT} python /tools/ldsc/ldsc.py --help'''
    out = subprocess.run(call.split(' '), capture_output=True)
    assert out.returncode == 0


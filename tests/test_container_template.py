# encoding: utf-8

"""
Test module for ``container_template.sif`` singularity build 
or ``container_template`` dockerfile build

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
    pth = os.path.join('containers', 'container_template.sif')
    out = subprocess.run('singularity')
    cwd = os.getcwd()
    PREFIX = f'singularity run {pth}'
    PREFIX_MOUNT = PREFIX_MOUNT = f'singularity run --home={cwd}:/home/ {pth}'
except FileNotFoundError:    
    try:
        out = subprocess.run('docker')
        pwd = os.getcwd()
        PREFIX = f'docker run -p {port}:{port} container_template'
        PREFIX_MOUNT = (f'docker run -p {port}:{port} ' + 
            f'--mount type=bind,source={pwd},target={pwd} container_template')
    except FileNotFoundError:
        raise FileNotFoundError('Neither `singularity` nor `docker` found in PATH. Can not run tests!')

def test_assert():
    """dummy test that should pass"""
    assert True

def test_container_template_python():
    """test that the Python installation works"""
    call = f'{PREFIX} python --version'
    out = subprocess.run(call.split(' '))
    assert out.returncode == 0

def test_container_template_python_script():
    pwd = os.getcwd() if PREFIX.rfind('docker') >= 0 else '.'
    call = f'''{PREFIX_MOUNT} python {pwd}/tests/extras/hello.py'''
    out = subprocess.run(call.split(' '), capture_output=True)
    assert out.returncode == 0
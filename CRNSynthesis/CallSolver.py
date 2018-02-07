import subprocess as sub
import re
import io
import os
from os import walk
import paramiko


def callLocalSolver(executableString):
    p = sub.Popen(executableString.split(), stdout=sub.PIPE, stderr=sub.PIPE)
    output, errors = p.communicate()
    return output


def callRemoteServer(executableString, server='', username='', password=''):
    p = sub.Popen(executableString.split(), stdout=sub.PIPE, stderr=sub.PIPE)
    ssh = paramiko.SSHClient()
    ssh.connect(server, username=username, password=password)
    ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command(cmd_to_execute)
    output, errors = p.communicate()
    return output


def ftpFileToServer():
    pass


def saveFileLocally():
    pass


import subprocess
import sys
import string
import os
import getpass
import time
from stat import *

################################################################################

PATH_SEQAN_GCC = '~/MyDocuments/VELOP/seqan2/'
PATH_SEQAN_WINDOWS = '../../'


################################################################################
### MODES = ['<mode> <target>', ...]
MODES = ['Simple justcompile']    
#MODES = ['Simple justcompile', 'Simple justrun']
#MODES = ['Release build', 'Release runonly']


### COMPILERS = [['<plaform>', '<name>', '<cmdline options>'], ...]
COMPILERS = [
#    ['gcc',     'g++-3.0',   'Compiler=g++-3.0'], 
#    ['gcc',     'g++-3.2',   'Compiler=g++-3.2'], 
    ['gcc',     'g++-3.3',   'Compiler=g++-3.3'], 
    ['gcc',     'g++-3.4',   'Compiler=g++-3.4'], 
    ['gcc',     'g++-4.1.2',   'Compiler=g++-4.1'], 
#    ['gcc',     'g++-4.1.1', 'Compiler=/import/testing/bin/g++'],
    ['windows', 'vc++-2003', 'Version=7'],
#    ['windows', 'vc++-2005', 'Version=8 "VSInstallDir=D:\\Program Files\\Microsoft Visual Studio 8\\ "']
]

### SERVERS = [['<plaform>', '<server>'], ...]
SERVERS = [
#    ['gcc', 'mouse'], 
#    ['gcc', 'frosch'], 
#    ['gcc', 'kiefer'], 
#    ['gcc', 'fichte'], 
#    ['gcc', 'aplysia'], 
    ['gcc', 'nawab'], 

#    ['gcc', 'duesseldorf'], 
#    ['gcc', 'hamm'], 
#    ['gcc', 'koeln'], 
#    ['gcc', 'bielefeld'], 
#    ['gcc', 'duisburg'], 
#    ['gcc', 'muenster'], 
#    ['gcc', 'aachen'], 
#    ['gcc', 'dortmund'], 
#    ['gcc', 'essen'], 
#    ['gcc', 'hagen'], 

    ['gcc', 'guangzhou'], 
    ['gcc', 'harbin'], 
    ['gcc', 'chongqing'], 
    ['gcc', 'shenyang'], 
    ['gcc', 'wuhan'], 
    ['gcc', 'chengdu'], 
    ['gcc', 'tianjin'], 
    ['gcc', 'xian'], 
    ['gcc', 'peking'], 

#    ['gcc', 'leningrad'], 
#    ['gcc', 'jekatarinburg'], 
#    ['gcc', 'irkutsk'], 
#    ['gcc', 'wladiwostok'], 
#    ['gcc', 'kaliningrad'], 
#    ['gcc', 'wolgograd'], 
#    ['gcc', 'omsk'], 
#    ['gcc', 'nowosibirsk'], 
#    ['gcc', 'moskau'], 

    ['windows', 'localhost']
]
    
LOAD_MAX = 1
SERVER_FAILURE_LIMIT = 3
 
################################################################################

PATH_MODULES = '../../projects/tests'
FILE_VIEWER = 'index.html'
FILE_OUTPUT = 'data.js'
FILE_OUTPUT_OLD = 'data_old.js'
PATH_RESULTFILES = 'results/'
PATH_GCC_RESULTFILES = PATH_SEQAN_GCC + 'misc/multitest/' + PATH_RESULTFILES

################################################################################

MODULES = []
ALL_STATES = ''
USERNAME = 'unknown'
PASSWORD = ''

################################################################################

def main():
    init()
    createPage()
    showPage()
    while not finishedAll():
        scheduleTasks()
        if pollTasks(): createPage()
        time.sleep(0.5)
    pollTasks()
    createPage()


################################################################################

def init():
    global SERVERS, COMPILERS, MODES, PATH_MODULES, LOAD_MAX
    global PATH_RESULTFILES
    global MODULES, TASKS, LOAD
    global USERNAME, PASSWORD
    
    USERNAME = getpass.getuser()
    PASSWORD = getpass.getpass()
    
    
    MODULES = []
    for f in os.listdir(PATH_MODULES):
        if (f == 'CVS') or (f == '.svn'): continue
        if f[0] == '_': continue
        p = PATH_MODULES + "/" + f
        m = os.stat(p)[ST_MODE]
        if S_ISDIR(m):
            MODULES.append(f)
            
    MODULES.sort()
    
    TASKS = []
    for mode in MODES:
        compnum = 0
        for [platform, compilername, compiler] in COMPILERS:
            compnum += 1
            for module in MODULES:
            	filename = module + '_' + platform + '_' + compilername + '_' + mode;
            	filename = filename.replace(' ', '_');
                task = {
                    'platform': platform,
                    'compilername': compilername,
                    'compiler': compiler,
                    'mode': mode,
                    'module': module,
                    'state': 'pending',
                    'filename': filename,
                    'proc' : None                    
                }
                TASKS.append(task)

    servers = []
    i = 0
    for [platform, servername] in SERVERS:
        server = {
            'i': i,
            'platform': platform,
            'server': servername,
            'load': 0,
            'loadmax': LOAD_MAX,
            'count': 0,
            'failures': 0
        }
        i += 1
        servers.append(server)
    SERVERS = servers
    
    if not os.access(PATH_RESULTFILES, os.F_OK):
        os.mkdir(PATH_RESULTFILES)

################################################################################
    
def javaScriptEscape(s):
    s = s.replace("'", "\\'")
    s = s.replace("\n", "'\n+ '")
    s = "DATA = '" + s + "';"
    return s;
    

def createPage():
    global SERVERS, COMPILERS, MODES, MODULES, TASKS
    global FILE_OUTPUT, FILE_OUTPUT_OLD

    s = ''
    i = 0
    for mode in MODES:
        s += '<div class=mode>' + mode + '</div>\n'
        s += '<table cellspacing=0 cellpadding=0 border=0>\n'
        s += '  <tr><th>&nbsp;</th>'
        for module in MODULES:
            s += '<th>' + module + '</th>'
        s += '  </tr>\n'
        
        for [platform, compilername, compiler] in COMPILERS:
            s += '  <tr>'
            s += '<td><nobr>' + platform + ': ' + compilername + '</nobr></td>'
            for module in MODULES:
                task = TASKS[i]
                i += 1
                state = task['state']
                if state != 'pending': server = task['server']['server']
                else: server = ''

                link = '<a href="' + PATH_RESULTFILES + task['filename'] + '.txt">'
                
                if state == 'pending': c = '<td class=pending align=middle>&nbsp;</td>'
                if state == 'started': c = '<td class=started align=middle>' + link + '<i>(' + server + ')</i></a></td>'
                if state == 'error': c = '<td class=error align=middle>' + link + 'error</a></td>'
                if state == 'ok': c = '<td class=ok align=middle>' + link + 'ok</a></td>'
                s += c
                
            s += '</tr>\n'
            
        s += '</table>\n'
    
    s += '<div class=mode>Server loads</div>\n'
    s += '<table cellspacing=0 cellpadding=0 border=0>\n'
    s += '<tr><td>&nbsp;</td><th>active</th><th>finished</th><th>failures</th></tr>\n'
    for server in SERVERS:
        if server['failures'] >= SERVER_FAILURE_LIMIT: cls = 'class=failed'
        else: cls = ''
        s += '<tr><th ' + cls + '>' + server['server'] + '</td><td align=middle ' + cls + '>' + str(server['load']) + '</td><td align=middle ' + cls + '>' + str(server['count']) + '</td><td align=middle ' + cls + '>' + str(server['failures']) + '</td></tr>'
    s += '</table>'
    
    if os.access(FILE_OUTPUT, os.F_OK):
	    f = open(FILE_OUTPUT)
	    t = f.read()
	    f.close()
    else:
    	t = ''
    	
    f = open(FILE_OUTPUT_OLD, "wb")
    f.write(t)
    f.close()

    f = open(FILE_OUTPUT, "wb")
    f.write(javaScriptEscape(s))
    f.close()

################################################################################

def showPage():
    global FILE_VIEWER
    c = [FILE_VIEWER]
    subprocess.Popen(c, shell=True)

################################################################################
    
def finishedAll():
    global TASKS
    
    for task in TASKS:
        if (task['state'] == 'pending') or (task['state'] == 'started'): return False
    return True

################################################################################

def scheduleTasks():
    global TASKS, SERVERS, SERVER_FAILURE_LIMIT
    
    for server in SERVERS:
        if (server['load'] < server['loadmax']) and (server['failures'] < SERVER_FAILURE_LIMIT):
            startTaskForServer(server)
            
        
################################################################################
       
def pollTasks():
    global TASKS
    global ALL_STATES
    
    old_all_states = ALL_STATES
    ALL_STATES = 'S'
    refresh = False

    for task in TASKS:
        ALL_STATES += task['state']
        
        if task['state'] != 'started': continue
        
        p = task['proc']
        if p == None: continue
        
        if (p.poll() != None): 
            onTaskStopped(task)
            refresh = True
        
    return refresh or (ALL_STATES != old_all_states)

################################################################################

def onTaskStopped(task):
    global ALL_STATES, PATH_RESULTFILES
    
    ALL_STATES = ''     # rebuild page
    
    server = task['server']
    server['load'] -= 1

    proc = task['proc']
    is_error = (proc.poll() != 0)
    
    if is_error:
        #read output file
        f = open(PATH_RESULTFILES + task['filename'] + '.txt', "rb")
        out = f.read()
        f.close()
        
        #test for rescheduling
        if (server['server'] != 'localhost') and (out.find('compile projects') == -1):
            #reschedule task
            server['failures'] += 1
            task['state'] = 'pending'
            task['proc'] = None
            print server['server'] + ' failed'
            return
            
    server['count'] += 1
    
    if is_error: task['state'] = 'error'
    else: task['state'] = 'ok'
    
    
    print task['filename'] + ":",
    if is_error: print "ERROR"
    else: print "OK"

################################################################################

def startTaskForServer(server):
    global TASKS, ALL_STATES
    
    for task in TASKS:
        if (task['state'] == 'pending') and (task['platform'] == server['platform']):
            server['load'] += 1
            task['state'] = 'started'
            task['server'] = server
            writeTaskHeader(task)
            task['proc'] = startTask(task, server)
            ALL_STATES = '' #rebuild page
            return
            
    #no more tasks for that platform
    server['loadmax'] = server['load']
    
    
def startTask(task, server):
    platform = task['platform']
    if platform == 'gcc': return startTaskGCC(task, server)
    if platform == 'windows': return startTaskWindows(task, server)
    
################################################################################
   
def writeTaskHeader(task):
    global PATH_RESULTFILES
    
    f = open(PATH_RESULTFILES + task['filename'] + '.txt', "wb")
    f.write('module:   ' + task['module'] + '\n')
    f.write('platform: ' + task['platform'] + '\n')
    f.write('compiler: ' + task['compilername'] + '\n')
    f.write('options:  ' + str(task['compiler']) + '\n')
    f.write('mode:     ' + task['mode'] + '\n')
    f.write('server:   ' + task['server']['server'] + '\n')
    f.write('result:   ' + task['state'] + '\n')
    f.write('------------------------------------------------------------------\n')
    f.close()

################################################################################

def startTaskGCC(task, server):
    global USERNAME, PASSWORD, PATH_SEQAN_GCC, PATH_GCC_RESULTFILES

    cmd = 'bash -c "'
    cmd += 'make -s -C' + PATH_SEQAN_GCC;
    cmd += ' Platform=' + task['platform']
    cmd += ' Project=' + task['module']
    cmd += ' Mode=' + task['mode']
    cmd += ' BuildFolder=' + task['platform'] + '/' + task['compilername']
    cmd += ' ' + task['compiler']
    cmd += ' 1>>' + PATH_GCC_RESULTFILES + task['filename'] + '.txt'
    cmd += ' 2>&1'
    cmd += '"'
    
    #print cmd
    
    c = ["plink.exe", "-batch", "-ssh", "-l", USERNAME, "-pw", PASSWORD, server['server'], cmd]
    return subprocess.Popen(c, shell=False)
        

def startTaskWindows(task, server):
    global PATH_SEQAN_WINDOWS


#    cmd = ["visual_studio\\make", 
#            "-s", 
#            "-C" + PATH_SEQAN_WINDOWS,
#            "Platform=" + task['platform'],
#            "Project=" + task['module'],
#            "Mode=" + task['mode'],
#            "BuildFolder=" + task['platform'] + "/" + task['compilername']
#         ] + task['compiler'] #+ ["1>>" + PATH_RESULTFILES + task['filename'] + ".txt", "2>&1"]
#    print cmd;  
#    return subprocess.Popen(cmd, shell=True)
    
    cmd = 'make -s -C' + PATH_SEQAN_WINDOWS
    cmd += ' Platform=' + task['platform']
    cmd += ' Project=' + task['module']
    cmd += ' Mode=' + task['mode']
    cmd += ' BuildFolder=' + task['platform'] + '/' + task['compilername']
    cmd += ' ' + task['compiler']
    cmd += ' 1>>' + PATH_RESULTFILES + task['filename'] + '.txt'
    cmd += ' 2>&1'
    cmd += '"'
    
#    print cmd
    
    c = [cmd]
    return subprocess.Popen(c, shell=True)


################################################################################


main()


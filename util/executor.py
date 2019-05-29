#coding:utf-8
import multiprocessing
import subprocess
import os
import time


class Loacl(object):
    """ a class for job run in local platform """
    def __init__(self):
        pass

    def submit(self, command, parameters):
        '''Submits a job
        '''
        pass

    def run_log(content, showTime=True):
        ##### Display text to defined stream, add timestamp
        logStream = sys.stdout
        if isinstance(content, list):
            modifText = " ".join(content)
        else:
            modifText = str(content)

        if showTime == True:
            modifText += "    " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        else:
            return modifText

    def runSysCommandSubProcess(command, showInLog=True, justTest=False):
        if showInLog == True:
            run_log(command, showTime=False)
            if justTest == True:
                return
        p = Popen(" ".join(command), shell=True, stderr=PIPE)
        return p.stderr.read().decode()

class SGE(object):
    ''' The SGE executor'''
    def __init__(self, mail_user=None):
        '''Constructor'''
        self.running = []
        self.queue = "normal"  # Default queue name is "normal"
        self.mem = 1000  # Request 1GB as a default
        self.out_dir = os.path.expanduser("~/tmp")
        self.cnt = 1
        self.project = "anopheles"
        self.mail_options = "a"
        self.mail_user = mail_user
        self.max_proc = 1000
        self.hosts = []
        self.cpus = 1

    def clean_done(self):
        '''Removes dead processes from the running list.'''
        ongoing = []
        status_file = "/tmp/farm-%d" % (os.getpid())
        os.system("qstat > %s 2>/dev/null" % status_file)
        f = open(status_file)
        f.readline()  # header
        f.readline()  # header
        for l in f:
            toks = list(filter(lambda x: x != "", l.rstrip().split(" ")))
            ongoing.append(int(toks[0]))
        os.remove(status_file)
        my_dels = []
        for r_idx, p in enumerate(self.running):
            if p not in ongoing:
                my_dels.append(r_idx)
        my_dels.reverse()
        for my_del in my_dels:
            del self.running[my_del]

    def wait(self, for_all=False, be_careful=60):
        '''Blocks according to some condition

           Args:
               for_all: Also waits if there is *ANY* job running (i.e.
                    block/barrier)

               be_careful: Wait X secs before starting. This is because
                       tasks take time to go into the pool.
        '''
        time.sleep(be_careful)
        self.clean_done()
        if for_all:
            while len(self.running) > 0:
                self.clean_done()
                time.sleep(1)

    def submit(self, command, parameters="", my_dir=os.getcwd()):
        '''Submits a job'''
        job_file = "/tmp/job-%d.%d" % (os.getpid(), self.cnt)
        w = open(job_file, "w")
        w.write("%s %s\n" % (command, parameters))
        w.close()

        if self.mail_user is not None:
            mail = "-m %s -M %s" % (self.mail_options, self.mail_user)
        else:
            mail = ""
        while len(self.running) > self.max_proc:
            self.wait(be_careful=5)
        hosts = ""
        if len(self.hosts) > 0:
            hosts = " -q "
        for host in self.hosts:
            hosts += "\*@%s" % host
            if host != self.hosts[-1]:
                hosts += ","
        job = "qsub %s %s -S /bin/bash -V -P %s -cwd -l h_vmem=%dm %s " % (
            mail, hosts, self.project, self.mem, job_file)
        status_file = "/tmp/farm-%d.%d" % (os.getpid(), self.cnt)
        os.system(job + " >%s 2>/dev/null" % status_file)
        f = open(status_file)
        l = f.readline()
        job = int(l.split(" ")[2])
        f.close()
        os.remove(status_file)
        os.remove(job_file)
        self.cnt += 1

        self.running.append(job)
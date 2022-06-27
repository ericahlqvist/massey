

import subprocess, threading

class Command(object):
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None

    def run(self, timeout):
        def target():
            
            self.process = subprocess.Popen(self.cmd, shell=True)
            self.process.communicate()
            #print('Thread finished')

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            #print('Terminating process')
            self.process.terminate()
            thread.join()
            return self.process.returncode


p = "3" # A prime > 2
Dmod = "8" # 3, 7, 4, 8
mod = ""
if (Dmod == "3" or Dmod == "7"):
    mod = "8"
else:
    mod = "16"

open_file = ""

if p == "3":
    open_file = "p_"+p+"_cyc_9_9_disc_"+Dmod+"_mod_"+mod+".txt"
else:
    open_file = "p_"+p+"_disc_"+Dmod+"_mod_"+mod+".txt"

file = open("discriminants/"+open_file)
lines = file.readlines()
for line in lines:
    my_str = ''.join(map(str, line))
    command = Command("./main "+p+" "+my_str)
    code = command.run(timeout=200)
    
    # if code == -15:
    #     res_file = open("output/"+p+"_"+Dmod+"mod"+mod+".txt", "a")
    #     print(my_str)
    #     res_file.write("{\"p\": \""+p+"\", \"D\": \""+my_str.strip()+"\", \"Z-rk\": \"-\", \"K-cyc\": \"-\", \"Lx-cyc\": \"-\", \"Ly-cyc\": \"-\", \"ZM\": \"-\"},\n")
    
file.close()

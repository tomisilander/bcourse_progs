import os

g_env = Environment(CCFLAGS= '-Wall')
#g_env['ENV']['LD_RUN_PATH'] = os.path.abspath('src/common') # not good

# SUPPORT FOR THE DEBUG FLAG

debug = ARGUMENTS.get('debug', 0)
if int(debug):
    g_env.Append(CCFLAGS = ' -g')
else:   
    g_env.Append(CCFLAGS = ' -O3')


# CHECKING DEPENDENCIES

conf = Configure(g_env)

if not conf.CheckHeader('Judy.h', language='C'):
    print 'Did not find Judy.h, exiting!'
    Exit(1)

if not conf.CheckLib('Judy', autoadd=0):
    print 'Did not find libJudy.a or Judy.lib, exiting!'
    Exit(1)

g_env = conf.Finish()


Export('g_env')

SConscript(['src/SConstruct'])

import sys
import GlobalVars, Diagnostics
import Library
import Fields
import Bulletin
import PlotField
import Functions
from Diagnostics import Launch

def spasso():
    ################# Welcome to SPASSO
    Library.welcome_message()
    Library.printMainMessage('Initialization')
    cruise = Library.choose_cruise()
    
    GlobalVars.init_dataDir(cruise)
    GlobalVars.init_prodDate()
    Library.Listproducts()

############################### GET DATA AND PLOT
    Library.printMainMessage("Get data and make figures")
    Library.tic()
    #get products and options
    products = [str(x) for x in GlobalVars.config.get('products','products').split(',')]
    opt = [str(x) for x in GlobalVars.config.get('plot_options','options').split(',')]
    mode = GlobalVars.config.get('cruises','mode')
    outmode = GlobalVars.config.get('cruises','outmode')
    
    for pr in products:
        nprod = GlobalVars.config.get('products',pr+'prod')
        #download data
        Library.printMessage("Downloading "+nprod)
        eval('Fields.'+nprod+'.download()')
        if mode == 'DT':
            date = GlobalVars.all_dates['datec_'+pr.lower()]
            if outmode == 'clim':
                Functions.climatology(nprod,date)
                
        #plot field
        Library.printMessage("Ploting "+nprod)
        PlotField.PlotField.Plot(cruise,pr,opt,type=None)
    Library.toc('satellite figures')
    
################################ DIAGNOSTICS
    Library.printMainMessage("Computing diagnostics")

    ### Eulerian
    Library.tic()
    if GlobalVars.Eul['diag'] is not None:
        if outmode == 'clim':
            Diagnostics.Launch(cruise,'eulerian',clim='on')
        else:
            Diagnostics.Launch(cruise,'eulerian')
    Library.toc('Eulerian diagnostics')
    ### Lagrangian
    if GlobalVars.Lag['diag'] is not None:
        if outmode == 'clim':
            Diagnostics.Launch(cruise,'lagrangian',clim='on')   
        else:
            Diagnostics.Launch(cruise,'lagrangian')
    
    ############################### COPYING FIGURE AND SAVED
    Library.printMainMessage("COPYING AND ZIPPING DATA")
    Library.cleantmp()
    Library.copyfiles()
    
    ############################### CREATE BULLETIN
    if GlobalVars.Bull['authors']!=None:
        Library.printMainMessage("CREATING BULLETIN")
        Bulletin.create(cruise)
    
    # ############################### SENDING EMAIL
    if GlobalVars.Email['sender']!=None:
        Library.printMainMessage("SENDING EMAIL")
        Library.send_email(cruise)

    return 0

############################### PROGRAM STARTS HERE
if __name__ == '__main__':

################# Usage

    args = sys.argv

    if (len(args) > 2):
        Library.usage()
        sys.exit(2)

################# Initilisation et récupération of 'config.ini' file

    cruise = args[1] # recupere le config ini de la cruise
    GlobalVars.configIni(cruise)
    Library.clean_wrk()
    GlobalVars.FigParam()
    GlobalVars.DiagParam()
    GlobalVars.EmailParam()
    GlobalVars.BulletinParam()
    GlobalVars.LibrariesPaths()
    
################# Get program starting time
    GlobalVars.init_date()

################# lauching the program with: spasso()
    if (spasso() == 0):
        Library.exit_program(0)
    else:
        Library.exit_program(1)

############################### end of program  ###############################

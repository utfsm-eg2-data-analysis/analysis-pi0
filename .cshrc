#ident  "@(#).cshrc     ver 1.0     Aug 20, 1996"
source /site/env/syscshrc
#source /group/clas/builds/PRODUCTION/packages/cms/jlab.cshrc DEVELOPMENT
source /group/clas/builds/PRODUCTION/packages/cms/jlab.cshrc PRODUCTION
setenv CLASTOOL /home/$LOGNAME/packages/ClasTool
setenv ANALYSIS /home/$LOGNAME/analyser
set history=100
set savehist=50
set noclobber
umask 077
use root
setenv LD_LIBRARY_PATH .:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH  ${CLASTOOL}/slib/${OS_NAME}:${LD_LIBRARY_PATH}
setenv LD_LIBRARY_PATH  ${ANALYSIS}/slib/:${LD_LIBRARY_PATH}
set prompt = "%m%b:/%C1[\!]> "
setenv LS_COLORS "*.cxx=31:*.C=31:*.cc=31:*.h=31:ex=32:di=34:ow=34:*rc=37:*akefile=1:*.root=6;1:*.mu=1;35:*.bos=36;1:*.gif=93:*.pdf=96;1:*.mu=35"
            alias       rm      rm -i
	                alias       cp      cp -i
			setenv SVNUSER https://jlabsvn.jlab.org/svnroot/clas/users/brooksw

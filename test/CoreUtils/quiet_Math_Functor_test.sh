#!/usr/bin/env bash
#
# Math/Functor.h uses std::auto_ptr, which triggers a warning
#
# Usage:  MathFunctor_wrapper_test.sh [include directories]
#
#

### BEGIN POSIX COMPLIANT ####
# Ubuntu 16 insists in using a non-bash default shell, and ctest insists to use the default shell;
# this horrible hack saves the day while issue #17234 is being addressed.
if [ -z "$BASH_VERSION" ] && [ "$1" != 'shellSwitch' ]
then
	echo "Attempting to switch to bash."
	bash -- "$0" 'shellSwitch' "$@"
	exit
fi
[ "$1" = 'shellSwitch' ] && shift # here bash would be ok...
### END POSIX COMPILANT ######


declare -a IncludeDirectives
for Dir in "$@" ; do
	IncludeDirectives=( "${IncludeDirectives[@]}" "-I${Dir}" )
done

Cmd=( g++ -std='c++14' -Wall -Werror=deprecated-declarations -pedantic "${IncludeDirectives[@]}" -x'c++' -c )
cat <<EOH
===============================================================================
CMD> ${Cmd[@]}
===============================================================================
EOH
"${Cmd[@]}" - <<EOC
#include "Math/Functor.h"
EOC
res=$?
cat <<EOM
-------------------------------------------------------------------------------
EXIT CODE: ${res} (errors are expected in the compilation above)
EOM
if [[ $res == 0 ]]; then
	echo "Unexpected success in compiling 'Math/Functor.h'. May be time for dropping its wrapper." >&2
	return 1
fi

# this should work instead:
cat <<EOH
===============================================================================
CMD> ${Cmd[@]}
===============================================================================
EOH
Cmd=( g++ -std='c++14' -Wall -Werror=deprecated-declarations -pedantic "${IncludeDirectives[@]}" -x'c++' -c )
"${Cmd[@]}" - <<EOC
#include "larcorealg/CoreUtils/quiet_Math_Functor.h"
EOC
res=$?
cat <<EOM
-------------------------------------------------------------------------------
EXIT CODE: ${res} (errors are NOT expected in the compilation above)
===============================================================================
EOM


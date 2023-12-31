#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | www.openfoam.com
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2011-2016 OpenFOAM Foundation
#     Copyright (C) 2017-2023 OpenCFD Ltd.
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, distributed under GPL-3.0-or-later.
#
#------------------------------------------------------------------------------
# Script
#     doxyFilter.sed
#
# Description
#     Transform human-readable tags such as 'Description' into the Doxygen
#     equivalent
#------------------------------------------------------------------------------

# New FSF address
/^License/,/\*\//{
/^License/,\%http://www.gnu.org/licenses%{
s?^License.*?\*\/\
\/\*! \\file %realFilePath%\
<b>Original source file</b> <a href="%filePath%">%fileName%</a>\
\
?
/^    /d
}


# Remove entry
/^Application *$/{
N
N
d
}


# Remove entry
/^Global *$/{
N
N
d
}


# Primitive
#     typename
# =>
# \\relates typename
#
/^Primitive *$/,/^[^ ]/{
s/^Primitive *$//
s/^    /\\relates /
}


# Class
#     Foam::className
# =>
# \\class Foam::className
#
# Class
#     Foam::namespaceName::
#         className
# =>
# \\class Foam::namespaceName::className
#
/^Class *$/{
N
:loop
/.*:: *$/{
N
s/^ *\(.*\) *\n *\(.*\) */\1\2/
}
t loop
s/Class *\n *\(.*\) */\\class \1/
}


# Group
#     groupName
# =>
# \ingroup groupName
#
/^Group *$/,/^[^ ]/{
s/^Group//
s/^    /\\ingroup /
}


# Namespace
#     namespaceName
# =>
# \namespace namespaceName
#
/^Namespace *$/,/^[^ ]/{
s/^Namespace//
s/^    /\\namespace /
}


# Typedef
#     Foam::def
# =>
# \typedef Foam::def
/^Typedef *$/,/^[^ ]/{
s/^Typedef//
s/^    /\\typedef /
}


# Add anchor and use \brief
# the first paragraph will be 'brief' and the others 'detail'
/^Description *$/,/^[^ ]/{
/^Description/c\
<a class="anchor" name="Description"></a>\
\\brief
s/^    //
}

/^Usage *$/,/^[^ ]/{
/^Usage/c\
\\par Usage
s/^    //
}


/^Environment *$/,/^[^ ]/{
/^Environment/c\
\\par Environment
s/^    //
}


/^See *[Aa]lso *$/,/^[^ ]/{
/^See *[Aa]lso/c\
\\see
s/^    //
}

/^Note *$/,/^[^ ]/{
/^Note/c\
\\note
s/^    //
}


/^Warning *$/,/^[^ ]/{
/^Warning/c\
\\warning
s/^    //
}


/^SourceFiles *$/,/^$/{
s?SourceFiles?\\par Source files\
<ul>\
  <li><a href="%filePath%">%fileName%</a></li> ?
s? *\([a-zA-Z0-9]*\.[a-zA-Z]*\)?  <li><a href="%dirName%/\1">\1</a></li>?
s?^$?</ul>?
}


/fileName%<\/a><\/li>$/{
N
s?\n$?</ul>?g
s/<\/li>\n/<\/li> /
s? *\([a-zA-Z0-9]*\.[a-zA-Z]*\)?  <li><a href="%dirName%/\1">\1</a></li>?
}

s/.*\*\//\*\//


# Convert \heading in source files to bold font and add some space
s#\\heading \(.*\)#<br><b>\1</b>#g

# Add a linebreak
s#\\linebreak#<br>#g

}


#------------------------------------------------------------------------------

#!/usr/bin/env ruby
require 'optparse'
################################################################
# CONFIGURATION:
CFTypeInfo             = {
  'int'                  => {:namedConst       => 'c_int'                , :ftype => 'integer'},
  'short int'            => {:namedConst       => 'c_short'              , :ftype => 'integer'},
  'long int'             => {:namedConst       => 'c_long'               , :ftype => 'integer'},
  'long long int'        => {:namedConst       => 'c_long_long'          , :ftype => 'integer'},
  'signed char'          => {:namedConst       => 'c_signed_char'        , :ftype => 'integer'},
  'unsigned char'        => {:namedConst       => 'c_signed_char'        , :ftype => 'integer'},
  'size_t'               => {:namedConst       => 'c_size_t'             , :ftype => 'integer'},
  'int8_t'               => {:namedConst       => 'c_int8_t'             , :ftype => 'integer'},
  'int16_t'              => {:namedConst       => 'c_int16_t'            , :ftype => 'integer'},
  'int32_t'              => {:namedConst       => 'c_int32_t'            , :ftype => 'integer'},
  'int64_t'              => {:namedConst       => 'c_int64_t'            , :ftype => 'integer'},
  'int_fast8_t'          => {:namedConst       => 'c_int_fast8_t'        , :ftype => 'integer'},
  'int_fast16_t'         => {:namedConst       => 'c_int_fast16_t'       , :ftype => 'integer'},
  'int_fast32_t'         => {:namedConst       => 'c_int_fast32_t'       , :ftype => 'integer'},
  'int_fast64_t'         => {:namedConst       => 'c_int_fast64_t'       , :ftype => 'integer'},
  'int_least8_t'         => {:namedConst       => 'c_int_least8_t'       , :ftype => 'integer'},
  'int_least16_t'        => {:namedConst       => 'c_int_least16_t'      , :ftype => 'integer'},
  'int_least32_t'        => {:namedConst       => 'c_int_least32_t'      , :ftype => 'integer'},
  'int_least64_t'        => {:namedConst       => 'c_int_least64_t'      , :ftype => 'integer'},
  'intmax_t'             => {:namedConst       => 'c_intmax_t'           , :ftype => 'integer'},
  'intptr_t'             => {:namedConst       => 'c_intptr_t'           , :ftype => 'integer'},

  'float'                => {:namedConst       => 'c_float'              , :ftype => 'real'},
  'double'               => {:namedConst       => 'c_double'             , :ftype => 'real'},
  'long double'          => {:namedConst       => 'c_long_double'        , :ftype => 'real'},

  'float _Complex'       => {:namedConst       => 'c_float_complex'      , :ftype => 'complex'},
  'double _Complex'      => {:namedConst       => 'c_double_complex'     , :ftype => 'complex'},
  'long double _Complex' => {:namedConst       => 'c_long_double_complex', :ftype => 'complex'},
  '_Bool'                => {:namedConst       => 'c_bool'               , :ftype => 'logical'},
  'char'                 => {:namedConst       => 'c_char'               , :ftype => 'character'}
}
# how the module should be invoked from fortran
ModuleName    = 'mo_cdi'
# which conversion is to use generating the fortran routine names, this could
# be any ruby String method.
FNameMethod   = :noop
# FNameMethod = :downcase
# FNameMethod = :noop
# FNameMethod = :upcase
class String;def noop;self;end;end
# fortran subroutines are not allowed to have a parameters with the same name,
# so in case of a match, these parameters have to get an new name
FParamExtension = 'v'
# Naming convention from above: 
# all non scalar variables should have the postfix '_vec' in there name
Vectors = /_vec$/i
################################################################################
FortranMaxLineLength = 132
# read the c header file and grep out name, return type and paramterlist of
# each function
def getFuncInfo(filename)
  typelist = %w[char int float double void]
  cppflags = ENV['CPPFLAGS'].nil? ? '' : ENV['CPPFLAGS']
  funclist = IO.popen("cpp #{cppflags} #{filename} | cpp -fpreprocessed").readlines.delete_if {|line| line.include?('#')}.collect {|line| line.chomp}
  # delete everything, that do not look like a function prototype
  typeRegexp = /^.*(#{typelist.join('|')}) \**\w+\s*\(.*\)/
  funclist.delete_if {|line|
    not typeRegexp.match(line.lstrip)
  }
  funclist.collect! {|line|
    md = /(\w+)+ +(\**)(\w+)\s*\((.*)\)/.match(line)
    returnType, returnPointer, funcName, paramList = md[1,4]
    paramList = paramList.split(',').collect {|p| p.split(' ').each {|_p| _p.strip}}
    [funcName, returnType, returnPointer, paramList]
  }
  funclist
end

# grep the #define C-Constants, which should be available within the fortran CDI API
def getDefines(filename)
  defines = File.open(filename,'r').readlines.grep(/^#define/).collect {|line|
    md = / +(\w+) +(-*\d+)/.match(line)
  }.select {|item| not item.nil?}.collect {|match| match[1..2]}
end

# create continuation for lines longer that 132 sign, which would create an
# error with some fortran compilers
def genContinuation(iline)
  # leave the input untouched
  line = iline
  # try to create readable line breaks, i.e. do not split name, labels or keywords
  regexp = /,[^,]*$/

  matchIndex = line[0,FortranMaxLineLength] =~ regexp

  if matchIndex.nil?
    line[FortranMaxLineLength-2] = "&\n"+line[FortranMaxLineLength-2,1]+"&"
  else
    line[matchIndex] = ",&\n" + ' '*7
  end
  line
end


# create fortran module variables
def genModParams(vars)
  vars.collect {|var, value|
    "      integer, parameter :: #{var} = #{value}"
  }.join("\n")
end

# return the fortran version of the c parameters (name, named constant within
# iso-c-bindig, fortran type)
def fortranParamsWithTypes(paramList)
  paramList.collect {|paramInfo|
    ctype, param = paramInfo[-2,2]
    ftype, nc    = CFTypeInfo[ctype][:ftype], CFTypeInfo[ctype][:namedConst]
    # fortran parameters can be configured out of the c parameters
    [param.send(FNameMethod),nc,ftype, paramInfo.include?('const')]
  }
end
def startMod(name)
  "
module #{name}
      use, intrinsic :: iso_c_binding

      implicit none

      private
  "
end
def endMod(name)
  "\nend module #{name}\n"
end
def isBadFunction(returnType, returnPointer, paramList)
  return true if (
    # external return type
    (returnType != 'void' and not CFTypeInfo.keys.include?(returnType)) or
    # pointer2pointer return type
    returnPointer.length > 1
  )
  paramList.each {|paramInfo|
    next if paramInfo == ['void']
    next if paramInfo[0] == 'void' and /^\*\w+$/.match(paramInfo[1])
    ctype, param = paramInfo[-2,2]
    # TJ: unnamed arguments shouldn't be matched at all but because
    # pointer * is parsed as part of the name those need to be rejected here
    return true if (
      CFTypeInfo[ctype].nil? or                        # derived data type
      param == '*' or                                  # unnamed pointer
      param[0,2] == '**' or                            # pointer2pointer
      (param[0,1] == '*' and /\w\[\d*\]/.match(param)) # array of pointers
    )
  }
  return false
end

def hasDimension(paramName, paramFType)
  return true if paramFType == 'character'
  return true if ( %w[integer real].include?(paramFType) and Vectors.match(paramName) )
  return false
end

# collect further information about the original c type of the params and
# in case of a match between param name and function name, the parameter name
# is changed
def setFortranParams(paramswithtypes,fFuncname)
  # special treatment of empty parameter list
  return [[],[]] if paramswithtypes == [[]]

  originalParameters = paramswithtypes.transpose[0]
  paramswithtypes.collect {|param, paramCType, fType, isConstant|
    # test for pointers/arrays
    isPointer, isArray, arraySize = false, false, nil
    if param[0,1] == '*'
      isPointer = true
      # remove '*' from funcnames and paramnames
      param.sub!('*','')
    end
    if ( md = /\w\[(\d*)\]/.match(param); not md.nil? )
      isArray   = true
      arraySize = md[1] == '' ? '*' : md[1]
      param.tr!("[#{md[1]}]",'')
    end

    # change param name if it equals the funcname
    if param == fFuncname or param.downcase == fFuncname.downcase
      param += FParamExtension
      # but maybe the result is the name of another parameter. -> play it again sam...
      param += FParamExtension while ( originalParameters.include?(param) )
    end
    if param[0,1] == '_'
      param = 'p' + param
    end

    [param,paramCType,fType,isPointer,isArray,arraySize,isConstant]
  }
end

def printParams(fParams, indent)
  out = ''
  fParams.each {|param,paramType,ftype,ispointer,isarray,arraysize,isconstant|
    dimension     = isarray ? "dimension(#{arraysize})" : ( (ispointer and hasDimension(param, ftype) ) ? 'dimension(*)' : nil)
    intent, value = nil, nil
    if (ispointer or isarray)
      if not isconstant
        intent = 'intent(out)'
      else
        intent = 'intent(in)'
      end unless paramType == 'c_char'
    else
      #intent = 'intent(in)'
      value  = 'value'
    end

    typeinfo = [value,intent,dimension].select {|s| ! s.nil?}.join(',')
    out << "  #{indent}#{ftype}"
    out << (paramType == 'c_ptr' ? '' : "(kind=#{paramType})")
    out << ", #{typeinfo} :: #{param}\n"
  }
  return out
end

# creates the actual binding within the module for the given c function
# unsupported types of function a ignored, see RESTRICTIONS (top) for details
def genInterface(cFuncname, returnType, returnPointer, paramList)

  # do not create interfaces for unsuppoerted function
  if isBadFunction( returnType, returnPointer, paramList)
    warn "parameterlist of '#{cFuncname}' is not supported -> function ignored."
    return ['','']
  end
  return ['', ''] if (cFuncname[0,1] == '_')
  # the void argument type can be left out: if 'void' occures in the
  # parameterlist (which is only the case, if it is the only parameter
  # information), it is simply removed and a empty paramter list is left.
  paramList = [[],[]] if paramList.flatten == [ 'void' ] or paramList.flatten.empty?

  out = ''
  isWrapper = false

  # create new names for the fortran routines, see CONFIGURATION (top)
  fFuncname = cFuncname.send(FNameMethod)

  paramsWithTypes = []
  # divide between empty and non empty parameter lists
  if paramList == [[],[]]
    fParams, fParams4Import, fTypes4Import = [], [],[]
  else
    # collect information for setting the correct fortran type for each parameter
    paramsWithTypes = paramList.collect {|paramInfo|
      ctype, param = paramInfo[-2,2]
      ptr_match = /^\*(\w+)$/.match(param)
      if (/^(:?const)? *void$/.match(ctype) and ptr_match)
        param = ptr_match[1]
        ftype = 'type(c_ptr)'
        nc = 'c_ptr'
      else
        ftype, nc    = CFTypeInfo[ctype][:ftype], CFTypeInfo[ctype][:namedConst]
      end
      # fortran parameters can be configured out of the c parameters
      [param.send(FNameMethod),nc,ftype, paramInfo.include?('const')]
    }

    # get deeper information about parameterlists and filter out functions with
    # non supported parameters
    fParams = setFortranParams(paramsWithTypes,fFuncname)
    fParams4Import, fTypes4Import, = fParams.transpose
  end

  fReturnInfo                    = CFTypeInfo[returnType]
  indent = '      '
  if not fReturnInfo.nil?
    fReturnString = "function"
    fEnd    = ["#{fReturnInfo[:ftype]}(kind=#{fReturnInfo[:namedConst]}) :: #{fFuncname}\n", "end function" ]
    fTypes4Import << fReturnInfo[:namedConst]
  elsif returnType == 'void'
    fReturnString = "subroutine"
    fEnd    = ["end subroutine"]
  end
  fFuncname_suffix = ''

  if (returnType == 'char' and returnPointer == '*')
    out << "#{indent}function #{fFuncname}(#{fParams4Import.join(',')})\n"
    out << printParams(fParams, indent)
    fTypes4Import.unshift('c_ptr') unless fTypes4Import.include?('c_ptr')
    isWrapper = true
    fwEnd = [ "end function" ]
    fFuncname_suffix = '_c'
    indent = '        '
    fEnd = [ "type(c_ptr) :: #{fFuncname}#{fFuncname_suffix}\n",
             "end function" ]
  end
  out << "#{indent}interface
#{indent}  #{fReturnString} #{fFuncname}#{fFuncname_suffix}(#{fParams4Import.join(',')}) bind(c,name='#{cFuncname}')\n"
  out << "#{indent}    import :: #{fTypes4Import.uniq.join(',')}\n" unless fTypes4Import.empty?

  out << printParams(fParams, indent + '  ')

  fEnd.each_with_index do |line, i|
    extra_indent = '    '
    if i == fEnd.length - 1
      extra_indent = '  '
    end
    out << indent + extra_indent + line
  end
  out << " #{fFuncname}#{fFuncname_suffix}\n#{indent}end interface\n"
  if (returnType == 'char' and returnPointer == '*')
    indent = '      '
    out << "#{indent}  character(len=1, kind=c_char), pointer :: #{fFuncname}(:)
#{indent}  type(c_ptr) :: cptr
#{indent}  integer :: slen(1)

#{indent}  cptr = #{fFuncname}#{fFuncname_suffix}("
    out << fParams4Import.join(",&
#{indent}    ") << ")\n"
    out << "#{indent}  #{fFuncname} => null()\n"
    out << "#{indent}  slen(1) = int(strlen(cptr))\n"
    out << "#{indent}  call c_f_pointer(cptr, #{fFuncname}, slen)\n"
    fwEnd.each_with_index do |line,i|
      extra_indent = '  '
      if i == fwEnd.length - 1
        extra_indent = ''
      end
      out << indent + extra_indent + line
    end
    out << " #{fFuncname}\n"
  end
  [out, makePublic(fFuncname), isWrapper]
end
def makePublic(*fFuncnameList)
  fFuncnameList.collect {|fname| 
    "      public :: #{fname.tr('*','')}"
  }.join("\n") << "\n"
end
def ctrim
  "
    subroutine ctrim(str)
    character(kind=c_char), intent(inout) :: str(:)
    character(kind=c_char) :: c
    integer :: i

    do i=1,size(str)
      c = str(i)
      if (c == c_null_char) then
        str(i:size(str)) = ' '
        exit
      end if
    end do

    end subroutine ctrim\n"
end

def clen
  "
    function c_len(s) result(i)
      character(kind=c_char), intent(in) :: s(:)
      integer :: i
      do i = 1, size(s)
        if (s(i) == c_null_char) then
          exit
        end if
      end do
      i = i - 1
    end function\n"
end

################################################################################
if __FILE__ == $0
require 'optparse'
require 'pp'

debug = false
OptionParser.new do |opts|
  opts.on("-d","--debug") {debug = true}
  opts.on_tail("--help","-h","-H","Display this help message.") do
    puts <<-'END'
#== Synopsis
# Create Fortran iso-c-bindings form a given c header file
# 
#== Usage
#   binGen.rb <headerFile> <fortranLibraryFile> [<modName>] [--help|--debug]  
#
# headerFile:
#   A general c header file: function prototypes and '#defines' with numerical
#   value will be taken for the fortran module construction. Furthermore there
#   are some restrictions to what is provided in fortran: Currently internal
#   datatypes and (arrays|pointers) to internal datatypes are supported, i.e.
#   no arrays of pointers, no pointer to pointers, no typedefs
#
# fortranLibraryFile:
#   file name for generated bindings
#
# modName:
#   This will be the name of the Fortran module, so it has to obey the fortran 
#   restriction for module names. default: mo_cdi
#  
#== Author
# Ralf Mueller, ralf.mueller@zmaw.de
#
#== RESTRICTIONS:
# ONLY SUPPORT FOR INTERNAL DATATYPES AND (ARRAYS|POINTERS) TO INTERNAL
# DATATYPES, I.E. NO ARRAYS OF POINTERS, NO POINTER TO POINTERS, NO TYPEDEFS
#
#=== Special: naming convention
# Pointers can have different sizes, which cannot be detetermined by parsing a
# function prototype. Therefor a convention according the parameter names can
# be used to take this decision precisely:
# * Pointers to numbers are expected to be scalars unless the name of the
#   parameter end with '_vec'
# * Pointers to char are allways referenced as vectors
#
#== LICENSE: 
# BSD License
#
#  Copyright (c) 2009-2012, Ralf Mueller (ralf.mueller@zmaw.de)
#  All rights reserved.
#  
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#  
#      * Redistributions of source code must retain the above copyright notice,
#        this list of conditions and the following disclaimer.
#  
#      * Redistributions in binary form must reproduce the above copyright
#        notice, this list of conditions and the following disclaimer in the
#        documentation and/or other materials provided with the distribution.
#  
#      * The names of its contributors may not be used to endorse or promote
#        products derived from this software without specific prior written
#        permission.
#  
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
END
  exit
  end
end.parse!

  if ARGV[1] == nil
    warn 'no outputile given'
    exit
  end

  outputString      = ''
  modname           = ARGV[2].nil? ? ModuleName : ARGV[2]

  cDefines          = getDefines(ARGV[0])
  pp cDefines if debug
  unless cDefines.empty?
    makeModVarsPublic = makePublic(*cDefines.transpose[0]) 
    moduleVariables   = genModParams(cDefines)
  end

  interfaces, makepublics, subroutines = '', '', "contains\n"
  indent = '    '

  funcdecls = [ [ 'strlen', 'size_t', '', [ [ 'void', '*s' ] ] ] ]
  funcdecls.concat(getFuncInfo(ARGV[0]))

  funcdecls.each {| funcName, returnType, returnPointer, paramList|
    pp [funcName, returnType, returnPointer, paramList] if debug
    interface, makepublic, isWrapper = genInterface(funcName,returnType, returnPointer, paramList)
    if isWrapper
      subroutines << interface
    else
      interfaces  << interface
    end
    makepublics << makepublic
  }

  # add a specialized trim for wierd c strings
  makepublics << makePublic('ctrim') << makePublic('c_len')
  subroutines << ctrim << clen

  File.open(ARGV[1],"w") {|f|
    [ startMod(modname),
      moduleVariables ||= '',
      interfaces,
      makepublics,
      makeModVarsPublic ||= '',
      subroutines,
      endMod(modname)
    ].join("\n").split("\n").each {|line| 
      # check the length of each line before writing to file
      if line.length > FortranMaxLineLength
        f << genContinuation(line) << "\n"
      else
        f << line << "\n"
      end
    }
  }

end

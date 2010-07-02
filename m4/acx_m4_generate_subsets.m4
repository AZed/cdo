dnl Helper function for recursion
m4_define([_ACX_M4_GENERATE_SUBSETS],dnl
[m4_ifval([$1],dnl
[_ACX_M4_GENERATE_SUBSETS(m4_cdr($1),m4_unquote(m4_cdr($@)),[$3])dnl
m4_ifval([$2],dnl
[_ACX_M4_GENERATE_SUBSETS(m4_cdr($1),dnl
m4_dquote(m4_unquote($2,m4_car($1))),[$3])],dnl
[_ACX_M4_GENERATE_SUBSETS(m4_cdr($1),m4_dquote(m4_car($1)),[$3])])dnl
],[m4_if([$2],[],,[,])m4_dquote(m4_join([$3],m4_unquote(m4_car(m4_shift($@)))))])])
dnl
dnl ACX_M4_GENERATE_SUBSETS(SET,[SEPARATOR])
dnl
dnl generates list of all subsets of SET, where SET is a
dnl comma-seperated list of elements like [[A],[B],[C]]
dnl
AC_DEFUN([ACX_M4_GENERATE_SUBSETS],dnl
[m4_dquote(_ACX_M4_GENERATE_SUBSETS([$1],[],[$2]))],dnl
)])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:

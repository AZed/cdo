#ifndef _ERROR_H
#define _ERROR_H

#define  _FATAL     1     /* Error flag: exit on error  */
#define  _VERBOSE   2     /* Error flag: report errors  */
#define  _DEBUG     4     /* Error flag: debug          */

extern int _ExitOnError;  /* If set to 1, exit on error (default 1)       */
extern int _Verbose;      /* If set to 1, errors are reported (default 1) */
extern int _Debug;        /* If set to 1, debuggig (default 0)            */

void SysError(const char *caller, const char *fmt, ...);
void    Error(const char *caller, const char *fmt, ...);
void  Warning(const char *caller, const char *fmt, ...);
void  Message(const char *caller, const char *fmt, ...);

#endif  /* _ERROR_H */

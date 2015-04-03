/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2009 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef _PSTREAM_H
#define _PSTREAM_H

#define  streamOpenWrite     pstreamOpenWrite
#define  streamOpenRead      pstreamOpenRead
#define  streamOpenAppend    pstreamOpenAppend
#define  streamClose         pstreamClose

#define  streamInqFiletype   pstreamInqFiletype

#define  streamInqVlist      pstreamInqVlist
#define  streamDefVlist      pstreamDefVlist

#define  streamDefTimestep   pstreamDefTimestep
#define  streamInqTimestep   pstreamInqTimestep

#define  streamDefRecord     pstreamDefRecord
#define  streamInqRecord     pstreamInqRecord

#define  streamWriteRecord   pstreamWriteRecord
#define  streamReadRecord    pstreamReadRecord
/*
#define  streamCopyRecord    pstreamCopyRecord
*/

int     pstreamOpenWrite(const char *streamname, int filetype);
int     pstreamOpenRead(const char *streamname);
int     pstreamOpenAppend(const char *streamname);
void    pstreamClose(int pstreamID);

int     pstreamInqFiletype(int pstreamID);

void    pstreamDefVlist(int pstreamID, int vlistID);
int     pstreamInqVlist(int pstreamID);

void    pstreamDefTimestep(int pstreamID, int tsID);
int     pstreamInqTimestep(int pstreamID, int tsID);

void    pstreamDefRecord(int pstreamID, int  varID, int  levelID);
int     pstreamInqRecord(int pstreamID, int *varID, int *levelID);

void    pstreamWriteRecord(int pstreamID, double *data, int nmiss);
void    pstreamReadRecord(int pstreamID, double *data, int *nmiss);
void    pstreamCopyRecord(int pstreamIDdest, int pstreamIDsrc);

#endif  /* _PSTREAM_H */

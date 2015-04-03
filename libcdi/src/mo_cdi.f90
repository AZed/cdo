
module mo_cdi
      use, intrinsic :: iso_c_binding

      implicit none

      private
  
      integer :: CDI_UNDEFID = -1
      integer :: CDI_GLOBAL = -1
      integer :: CDI_BIGENDIAN = 0
      integer :: CDI_LITTLEENDIAN = 1
      integer :: CDI_REAL = 1
      integer :: CDI_COMP = 2
      integer :: CDI_BOTH = 3
      integer :: CDI_ESYSTEM = -10
      integer :: CDI_EINVAL = -20
      integer :: CDI_EUFTYPE = -21
      integer :: CDI_ELIBNAVAIL = -22
      integer :: CDI_EUFSTRUCT = -23
      integer :: CDI_EUNC4 = -24
      integer :: CDI_ELIMIT = -99
      integer :: FILETYPE_GRB = 1
      integer :: FILETYPE_GRB2 = 2
      integer :: FILETYPE_NC = 3
      integer :: FILETYPE_NC2 = 4
      integer :: FILETYPE_NC4 = 5
      integer :: FILETYPE_SRV = 6
      integer :: FILETYPE_EXT = 7
      integer :: FILETYPE_IEG = 8
      integer :: COMPRESS_NONE = 0
      integer :: COMPRESS_SZIP = 1
      integer :: COMPRESS_GZIP = 2
      integer :: COMPRESS_BZIP2 = 3
      integer :: COMPRESS_ZIP = 4
      integer :: COMPRESS_JPEG = 5
      integer :: DATATYPE_PACK = 0
      integer :: DATATYPE_PACK1 = 1
      integer :: DATATYPE_PACK2 = 2
      integer :: DATATYPE_PACK3 = 3
      integer :: DATATYPE_PACK4 = 4
      integer :: DATATYPE_PACK5 = 5
      integer :: DATATYPE_PACK6 = 6
      integer :: DATATYPE_PACK7 = 7
      integer :: DATATYPE_PACK8 = 8
      integer :: DATATYPE_PACK9 = 9
      integer :: DATATYPE_PACK10 = 10
      integer :: DATATYPE_PACK11 = 11
      integer :: DATATYPE_PACK12 = 12
      integer :: DATATYPE_PACK13 = 13
      integer :: DATATYPE_PACK14 = 14
      integer :: DATATYPE_PACK15 = 15
      integer :: DATATYPE_PACK16 = 16
      integer :: DATATYPE_PACK17 = 17
      integer :: DATATYPE_PACK18 = 18
      integer :: DATATYPE_PACK19 = 19
      integer :: DATATYPE_PACK20 = 20
      integer :: DATATYPE_PACK21 = 21
      integer :: DATATYPE_PACK22 = 22
      integer :: DATATYPE_PACK23 = 23
      integer :: DATATYPE_PACK24 = 24
      integer :: DATATYPE_PACK25 = 25
      integer :: DATATYPE_PACK26 = 26
      integer :: DATATYPE_PACK27 = 27
      integer :: DATATYPE_PACK28 = 28
      integer :: DATATYPE_PACK29 = 29
      integer :: DATATYPE_PACK30 = 30
      integer :: DATATYPE_PACK31 = 31
      integer :: DATATYPE_PACK32 = 32
      integer :: DATATYPE_CPX32 = 64
      integer :: DATATYPE_CPX64 = 128
      integer :: DATATYPE_FLT32 = 132
      integer :: DATATYPE_FLT64 = 164
      integer :: DATATYPE_INT8 = 208
      integer :: DATATYPE_INT16 = 216
      integer :: DATATYPE_INT32 = 232
      integer :: DATATYPE_UINT8 = 308
      integer :: DATATYPE_UINT16 = 316
      integer :: DATATYPE_UINT32 = 332
      integer :: DATATYPE_INT = 251
      integer :: DATATYPE_FLT = 252
      integer :: DATATYPE_TXT = 253
      integer :: DATATYPE_CPX = 254
      integer :: GRID_GENERIC = 1
      integer :: GRID_GAUSSIAN = 2
      integer :: GRID_GAUSSIAN_REDUCED = 3
      integer :: GRID_LONLAT = 4
      integer :: GRID_SPECTRAL = 5
      integer :: GRID_FOURIER = 6
      integer :: GRID_GME = 7
      integer :: GRID_TRAJECTORY = 8
      integer :: GRID_UNSTRUCTURED = 9
      integer :: GRID_CURVILINEAR = 10
      integer :: GRID_LCC = 11
      integer :: GRID_LCC2 = 12
      integer :: GRID_LAEA = 13
      integer :: GRID_SINUSOIDAL = 14
      integer :: ZAXIS_SURFACE = 0
      integer :: ZAXIS_GENERIC = 1
      integer :: ZAXIS_HYBRID = 2
      integer :: ZAXIS_HYBRID_HALF = 3
      integer :: ZAXIS_PRESSURE = 4
      integer :: ZAXIS_HEIGHT = 5
      integer :: ZAXIS_DEPTH_BELOW_SEA = 6
      integer :: ZAXIS_DEPTH_BELOW_LAND = 7
      integer :: ZAXIS_ISENTROPIC = 8
      integer :: ZAXIS_TRAJECTORY = 9
      integer :: ZAXIS_ALTITUDE = 10
      integer :: ZAXIS_SIGMA = 11
      integer :: ZAXIS_MEANSEA = 12
      integer :: TAXIS_ABSOLUTE = 1
      integer :: TAXIS_RELATIVE = 2
      integer :: TIME_CONSTANT = 1
      integer :: TIME_VARIABLE = 2
      integer :: TUNIT_SECOND = 1
      integer :: TUNIT_MINUTE = 2
      integer :: TUNIT_HOUR = 3
      integer :: TUNIT_DAY = 4
      integer :: TUNIT_MONTH = 5
      integer :: TUNIT_YEAR = 6
      integer :: TUNIT_QUARTER = 7
      integer :: TUNIT_3HOURS = 8
      integer :: TUNIT_6HOURS = 9
      integer :: TUNIT_12HOURS = 10
      integer :: TSTEP_INSTANT = 1
      integer :: TSTEP_AVG = 2
      integer :: TSTEP_ACCUM = 3
      integer :: TSTEP_MAX = 4
      integer :: TSTEP_MIN = 5
      integer :: TSTEP_DIFF = 6
      integer :: TSTEP_RANGE = 7
      integer :: TSTEP_INSTANT2 = 8
      integer :: TSTEP_INSTANT3 = 9
      integer :: CALENDAR_STANDARD = 0
      integer :: CALENDAR_PROLEPTIC = 1
      integer :: CALENDAR_360DAYS = 2
      integer :: CALENDAR_365DAYS = 3
      integer :: CALENDAR_366DAYS = 4
      integer :: CALENDAR_NONE = 5

      interface
        character(c_char) function cdiStringError(cdiErrno) bind(c,name='cdiStringError')
          import :: c_int,c_char
          integer(c_int), value :: cdiErrno
       end function cdiStringError
      end interface
  
      interface
        subroutine cdiDebug(debug) bind(c,name='cdiDebug')
          import :: c_int
          integer(c_int), value :: debug
       end subroutine cdiDebug
      end interface
  
      interface
        character(c_char) function cdiLibraryVersion() bind(c,name='cdiLibraryVersion')
          import :: c_char
       end function cdiLibraryVersion
      end interface
  
      interface
        subroutine cdiPrintVersion() bind(c,name='cdiPrintVersion')
       end subroutine cdiPrintVersion
      end interface
  
      interface
        subroutine cdiDefMissval(missval) bind(c,name='cdiDefMissval')
          import :: c_double
          real(c_double), value :: missval
       end subroutine cdiDefMissval
      end interface
  
      interface
        real(c_double) function cdiInqMissval() bind(c,name='cdiInqMissval')
          import :: c_double
       end function cdiInqMissval
      end interface
  
      interface
        subroutine cdiDefGlobal(string,val) bind(c,name='cdiDefGlobal')
          import :: c_char,c_int
          character(c_char), dimension(*) :: string
          integer(c_int), value :: val
       end subroutine cdiDefGlobal
      end interface
  
      interface
        subroutine cdiParamToString(param,paramstr,maxlen) bind(c,name='cdiParamToString')
          import :: c_int,c_char
          integer(c_int), value :: param
          character(c_char), dimension(*) :: paramstr
          integer(c_int), value :: maxlen
       end subroutine cdiParamToString
      end interface
  
      interface
        subroutine cdiDecodeParam(param,pnum,pcat,pdis) bind(c,name='cdiDecodeParam')
          import :: c_int
          integer(c_int), value :: param
          integer(c_int), intent(out) :: pnum
          integer(c_int), intent(out) :: pcat
          integer(c_int), intent(out) :: pdis
       end subroutine cdiDecodeParam
      end interface
  
      interface
        integer(c_int) function cdiEncodeParam(pnum,pcat,pdis) bind(c,name='cdiEncodeParam')
          import :: c_int
          integer(c_int), value :: pnum
          integer(c_int), value :: pcat
          integer(c_int), value :: pdis
       end function cdiEncodeParam
      end interface
  
      interface
        subroutine cdiDecodeDate(date,year,month,day) bind(c,name='cdiDecodeDate')
          import :: c_int
          integer(c_int), value :: date
          integer(c_int), intent(out) :: year
          integer(c_int), intent(out) :: month
          integer(c_int), intent(out) :: day
       end subroutine cdiDecodeDate
      end interface
  
      interface
        integer(c_int) function cdiEncodeDate(year,month,day) bind(c,name='cdiEncodeDate')
          import :: c_int
          integer(c_int), value :: year
          integer(c_int), value :: month
          integer(c_int), value :: day
       end function cdiEncodeDate
      end interface
  
      interface
        subroutine cdiDecodeTime(time,hour,minute,second) bind(c,name='cdiDecodeTime')
          import :: c_int
          integer(c_int), value :: time
          integer(c_int), intent(out) :: hour
          integer(c_int), intent(out) :: minute
          integer(c_int), intent(out) :: second
       end subroutine cdiDecodeTime
      end interface
  
      interface
        integer(c_int) function cdiEncodeTime(hour,minute,second) bind(c,name='cdiEncodeTime')
          import :: c_int
          integer(c_int), value :: hour
          integer(c_int), value :: minute
          integer(c_int), value :: second
       end function cdiEncodeTime
      end interface
  
      interface
        integer(c_int) function streamOpenRead(path) bind(c,name='streamOpenRead')
          import :: c_char,c_int
          character(c_char), dimension(*) :: path
       end function streamOpenRead
      end interface
  
      interface
        integer(c_int) function streamOpenWrite(path,filetype) bind(c,name='streamOpenWrite')
          import :: c_char,c_int
          character(c_char), dimension(*) :: path
          integer(c_int), value :: filetype
       end function streamOpenWrite
      end interface
  
      interface
        integer(c_int) function streamOpenAppend(path) bind(c,name='streamOpenAppend')
          import :: c_char,c_int
          character(c_char), dimension(*) :: path
       end function streamOpenAppend
      end interface
  
      interface
        subroutine streamClose(streamID) bind(c,name='streamClose')
          import :: c_int
          integer(c_int), value :: streamID
       end subroutine streamClose
      end interface
  
      interface
        subroutine streamSync(streamID) bind(c,name='streamSync')
          import :: c_int
          integer(c_int), value :: streamID
       end subroutine streamSync
      end interface
  
      interface
        subroutine streamDefVlist(streamID,vlistID) bind(c,name='streamDefVlist')
          import :: c_int
          integer(c_int), value :: streamID
          integer(c_int), value :: vlistID
       end subroutine streamDefVlist
      end interface
  
      interface
        integer(c_int) function streamInqVlist(streamID) bind(c,name='streamInqVlist')
          import :: c_int
          integer(c_int), value :: streamID
       end function streamInqVlist
      end interface
  
      interface
        integer(c_int) function streamInqFiletype(streamID) bind(c,name='streamInqFiletype')
          import :: c_int
          integer(c_int), value :: streamID
       end function streamInqFiletype
      end interface
  
      interface
        subroutine streamDefByteorder(streamID,byteorder) bind(c,name='streamDefByteorder')
          import :: c_int
          integer(c_int), value :: streamID
          integer(c_int), value :: byteorder
       end subroutine streamDefByteorder
      end interface
  
      interface
        integer(c_int) function streamInqByteorder(streamID) bind(c,name='streamInqByteorder')
          import :: c_int
          integer(c_int), value :: streamID
       end function streamInqByteorder
      end interface
  
      interface
        subroutine streamDefZtype(streamID,ztype) bind(c,name='streamDefZtype')
          import :: c_int
          integer(c_int), value :: streamID
          integer(c_int), value :: ztype
       end subroutine streamDefZtype
      end interface
  
      interface
        subroutine streamDefZlevel(streamID,zlevel) bind(c,name='streamDefZlevel')
          import :: c_int
          integer(c_int), value :: streamID
          integer(c_int), value :: zlevel
       end subroutine streamDefZlevel
      end interface
  
      interface
        integer(c_int) function streamInqZtype(streamID) bind(c,name='streamInqZtype')
          import :: c_int
          integer(c_int), value :: streamID
       end function streamInqZtype
      end interface
  
      interface
        integer(c_int) function streamInqZlevel(streamID) bind(c,name='streamInqZlevel')
          import :: c_int
          integer(c_int), value :: streamID
       end function streamInqZlevel
      end interface
  
      interface
        integer(c_int) function streamDefTimestep(streamID,tsID) bind(c,name='streamDefTimestep')
          import :: c_int
          integer(c_int), value :: streamID
          integer(c_int), value :: tsID
       end function streamDefTimestep
      end interface
  
      interface
        integer(c_int) function streamInqTimestep(streamID,tsID) bind(c,name='streamInqTimestep')
          import :: c_int
          integer(c_int), value :: streamID
          integer(c_int), value :: tsID
       end function streamInqTimestep
      end interface
  
      interface
        character(c_char) function streamFilename(streamID) bind(c,name='streamFilename')
          import :: c_int,c_char
          integer(c_int), value :: streamID
       end function streamFilename
      end interface
  
      interface
        character(c_char) function streamFilesuffix(filetype) bind(c,name='streamFilesuffix')
          import :: c_int,c_char
          integer(c_int), value :: filetype
       end function streamFilesuffix
      end interface
  
      interface
        integer(c_int) function streamNtsteps(streamID) bind(c,name='streamNtsteps')
          import :: c_int
          integer(c_int), value :: streamID
       end function streamNtsteps
      end interface
  
      interface
        subroutine streamReadVar(streamID,varID,data_vec,nmiss) bind(c,name='streamReadVar')
          import :: c_int,c_double
          integer(c_int), value :: streamID
          integer(c_int), value :: varID
          real(c_double), intent(out),dimension(*) :: data_vec
          integer(c_int), intent(out) :: nmiss
       end subroutine streamReadVar
      end interface
  
      interface
        subroutine streamWriteVar(streamID,varID,data_vec,nmiss) bind(c,name='streamWriteVar')
          import :: c_int,c_double
          integer(c_int), value :: streamID
          integer(c_int), value :: varID
          real(c_double), intent(in),dimension(*) :: data_vec
          integer(c_int), value :: nmiss
       end subroutine streamWriteVar
      end interface
  
      interface
        subroutine streamReadVarSlice(streamID,varID,levelID,data_vec,nmiss) bind(c,name='streamReadVarSlice')
          import :: c_int,c_double
          integer(c_int), value :: streamID
          integer(c_int), value :: varID
          integer(c_int), value :: levelID
          real(c_double), intent(out),dimension(*) :: data_vec
          integer(c_int), intent(out) :: nmiss
       end subroutine streamReadVarSlice
      end interface
  
      interface
        subroutine streamWriteVarSlice(streamID,varID,levelID,data_vec,nmiss) bind(c,name='streamWriteVarSlice')
          import :: c_int,c_double
          integer(c_int), value :: streamID
          integer(c_int), value :: varID
          integer(c_int), value :: levelID
          real(c_double), intent(in),dimension(*) :: data_vec
          integer(c_int), value :: nmiss
       end subroutine streamWriteVarSlice
      end interface
  
      interface
        subroutine streamInqRecord(streamID,varID,levelID) bind(c,name='streamInqRecord')
          import :: c_int
          integer(c_int), value :: streamID
          integer(c_int), intent(out) :: varID
          integer(c_int), intent(out) :: levelID
       end subroutine streamInqRecord
      end interface
  
      interface
        subroutine streamDefRecord(streamID,varID,levelID) bind(c,name='streamDefRecord')
          import :: c_int
          integer(c_int), value :: streamID
          integer(c_int), value :: varID
          integer(c_int), value :: levelID
       end subroutine streamDefRecord
      end interface
  
      interface
        subroutine streamReadRecord(streamID,data_vec,nmiss) bind(c,name='streamReadRecord')
          import :: c_int,c_double
          integer(c_int), value :: streamID
          real(c_double), intent(out),dimension(*) :: data_vec
          integer(c_int), intent(out) :: nmiss
       end subroutine streamReadRecord
      end interface
  
      interface
        subroutine streamWriteRecord(streamID,data_vec,nmiss) bind(c,name='streamWriteRecord')
          import :: c_int,c_double
          integer(c_int), value :: streamID
          real(c_double), intent(in),dimension(*) :: data_vec
          integer(c_int), value :: nmiss
       end subroutine streamWriteRecord
      end interface
  
      interface
        subroutine streamCopyRecord(streamIDdest,streamIDsrc) bind(c,name='streamCopyRecord')
          import :: c_int
          integer(c_int), value :: streamIDdest
          integer(c_int), value :: streamIDsrc
       end subroutine streamCopyRecord
      end interface
  
      interface
        subroutine streamInqGinfo(streamID,intnum,fltnum) bind(c,name='streamInqGinfo')
          import :: c_int,c_float
          integer(c_int), value :: streamID
          integer(c_int), intent(out) :: intnum
          real(c_float), intent(out) :: fltnum
       end subroutine streamInqGinfo
      end interface
  
      interface
        integer(c_int) function vlistCreate() bind(c,name='vlistCreate')
          import :: c_int
       end function vlistCreate
      end interface
  
      interface
        subroutine vlistDestroy(vlistID) bind(c,name='vlistDestroy')
          import :: c_int
          integer(c_int), value :: vlistID
       end subroutine vlistDestroy
      end interface
  
      interface
        integer(c_int) function vlistDuplicate(vlistID) bind(c,name='vlistDuplicate')
          import :: c_int
          integer(c_int), value :: vlistID
       end function vlistDuplicate
      end interface
  
      interface
        subroutine vlistCopy(vlistID2,vlistID1) bind(c,name='vlistCopy')
          import :: c_int
          integer(c_int), value :: vlistID2
          integer(c_int), value :: vlistID1
       end subroutine vlistCopy
      end interface
  
      interface
        subroutine vlistCopyFlag(vlistID2,vlistID1) bind(c,name='vlistCopyFlag')
          import :: c_int
          integer(c_int), value :: vlistID2
          integer(c_int), value :: vlistID1
       end subroutine vlistCopyFlag
      end interface
  
      interface
        subroutine vlistClearFlag(vlistID) bind(c,name='vlistClearFlag')
          import :: c_int
          integer(c_int), value :: vlistID
       end subroutine vlistClearFlag
      end interface
  
      interface
        subroutine vlistCat(vlistID2,vlistID1) bind(c,name='vlistCat')
          import :: c_int
          integer(c_int), value :: vlistID2
          integer(c_int), value :: vlistID1
       end subroutine vlistCat
      end interface
  
      interface
        subroutine vlistMerge(vlistID2,vlistID1) bind(c,name='vlistMerge')
          import :: c_int
          integer(c_int), value :: vlistID2
          integer(c_int), value :: vlistID1
       end subroutine vlistMerge
      end interface
  
      interface
        subroutine vlistPrint(vlistID) bind(c,name='vlistPrint')
          import :: c_int
          integer(c_int), value :: vlistID
       end subroutine vlistPrint
      end interface
  
      interface
        integer(c_int) function vlistNumber(vlistID) bind(c,name='vlistNumber')
          import :: c_int
          integer(c_int), value :: vlistID
       end function vlistNumber
      end interface
  
      interface
        integer(c_int) function vlistNvars(vlistID) bind(c,name='vlistNvars')
          import :: c_int
          integer(c_int), value :: vlistID
       end function vlistNvars
      end interface
  
      interface
        integer(c_int) function vlistNgrids(vlistID) bind(c,name='vlistNgrids')
          import :: c_int
          integer(c_int), value :: vlistID
       end function vlistNgrids
      end interface
  
      interface
        integer(c_int) function vlistNzaxis(vlistID) bind(c,name='vlistNzaxis')
          import :: c_int
          integer(c_int), value :: vlistID
       end function vlistNzaxis
      end interface
  
      interface
        subroutine vlistDefNtsteps(vlistID,nts) bind(c,name='vlistDefNtsteps')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: nts
       end subroutine vlistDefNtsteps
      end interface
  
      interface
        integer(c_int) function vlistNtsteps(vlistID) bind(c,name='vlistNtsteps')
          import :: c_int
          integer(c_int), value :: vlistID
       end function vlistNtsteps
      end interface
  
      interface
        integer(c_int) function vlistGridsizeMax(vlistID) bind(c,name='vlistGridsizeMax')
          import :: c_int
          integer(c_int), value :: vlistID
       end function vlistGridsizeMax
      end interface
  
      interface
        integer(c_int) function vlistGrid(vlistID,index) bind(c,name='vlistGrid')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: index
       end function vlistGrid
      end interface
  
      interface
        integer(c_int) function vlistGridIndex(vlistID,gridID) bind(c,name='vlistGridIndex')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: gridID
       end function vlistGridIndex
      end interface
  
      interface
        subroutine vlistChangeGridIndex(vlistID,index,gridID) bind(c,name='vlistChangeGridIndex')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: index
          integer(c_int), value :: gridID
       end subroutine vlistChangeGridIndex
      end interface
  
      interface
        subroutine vlistChangeGrid(vlistID,gridID1,gridID2) bind(c,name='vlistChangeGrid')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: gridID1
          integer(c_int), value :: gridID2
       end subroutine vlistChangeGrid
      end interface
  
      interface
        integer(c_int) function vlistZaxis(vlistID,index) bind(c,name='vlistZaxis')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: index
       end function vlistZaxis
      end interface
  
      interface
        integer(c_int) function vlistZaxisIndex(vlistID,zaxisID) bind(c,name='vlistZaxisIndex')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: zaxisID
       end function vlistZaxisIndex
      end interface
  
      interface
        subroutine vlistChangeZaxisIndex(vlistID,index,zaxisID) bind(c,name='vlistChangeZaxisIndex')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: index
          integer(c_int), value :: zaxisID
       end subroutine vlistChangeZaxisIndex
      end interface
  
      interface
        subroutine vlistChangeZaxis(vlistID,zaxisID1,zaxisID2) bind(c,name='vlistChangeZaxis')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: zaxisID1
          integer(c_int), value :: zaxisID2
       end subroutine vlistChangeZaxis
      end interface
  
      interface
        integer(c_int) function vlistNrecs(vlistID) bind(c,name='vlistNrecs')
          import :: c_int
          integer(c_int), value :: vlistID
       end function vlistNrecs
      end interface
  
      interface
        subroutine vlistDefTaxis(vlistID,taxisID) bind(c,name='vlistDefTaxis')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: taxisID
       end subroutine vlistDefTaxis
      end interface
  
      interface
        integer(c_int) function vlistInqTaxis(vlistID) bind(c,name='vlistInqTaxis')
          import :: c_int
          integer(c_int), value :: vlistID
       end function vlistInqTaxis
      end interface
  
      interface
        subroutine vlistDefTable(vlistID,tableID) bind(c,name='vlistDefTable')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: tableID
       end subroutine vlistDefTable
      end interface
  
      interface
        integer(c_int) function vlistInqTable(vlistID) bind(c,name='vlistInqTable')
          import :: c_int
          integer(c_int), value :: vlistID
       end function vlistInqTable
      end interface
  
      interface
        subroutine vlistDefInstitut(vlistID,instID) bind(c,name='vlistDefInstitut')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: instID
       end subroutine vlistDefInstitut
      end interface
  
      interface
        integer(c_int) function vlistInqInstitut(vlistID) bind(c,name='vlistInqInstitut')
          import :: c_int
          integer(c_int), value :: vlistID
       end function vlistInqInstitut
      end interface
  
      interface
        subroutine vlistDefModel(vlistID,modelID) bind(c,name='vlistDefModel')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: modelID
       end subroutine vlistDefModel
      end interface
  
      interface
        integer(c_int) function vlistInqModel(vlistID) bind(c,name='vlistInqModel')
          import :: c_int
          integer(c_int), value :: vlistID
       end function vlistInqModel
      end interface
  
      interface
        integer(c_int) function vlistDefVar(vlistID,gridID,zaxisID,timeID) bind(c,name='vlistDefVar')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: gridID
          integer(c_int), value :: zaxisID
          integer(c_int), value :: timeID
       end function vlistDefVar
      end interface
  
      interface
        subroutine vlistChangeVarGrid(vlistID,varID,gridID) bind(c,name='vlistChangeVarGrid')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: gridID
       end subroutine vlistChangeVarGrid
      end interface
  
      interface
        subroutine vlistChangeVarZaxis(vlistID,varID,zaxisID) bind(c,name='vlistChangeVarZaxis')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: zaxisID
       end subroutine vlistChangeVarZaxis
      end interface
  
      interface
        subroutine vlistInqVar(vlistID,varID,gridID,zaxisID,timeID) bind(c,name='vlistInqVar')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), intent(out) :: gridID
          integer(c_int), intent(out) :: zaxisID
          integer(c_int), intent(out) :: timeID
       end subroutine vlistInqVar
      end interface
  
      interface
        integer(c_int) function vlistInqVarGrid(vlistID,varID) bind(c,name='vlistInqVarGrid')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarGrid
      end interface
  
      interface
        integer(c_int) function vlistInqVarZaxis(vlistID,varID) bind(c,name='vlistInqVarZaxis')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarZaxis
      end interface
  
      interface
        integer(c_int) function vlistInqVarTime(vlistID,varID) bind(c,name='vlistInqVarTime')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarTime
      end interface
  
      interface
        subroutine vlistDefVarZtype(vlistID,varID,ztype) bind(c,name='vlistDefVarZtype')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: ztype
       end subroutine vlistDefVarZtype
      end interface
  
      interface
        integer(c_int) function vlistInqVarZtype(vlistID,varID) bind(c,name='vlistInqVarZtype')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarZtype
      end interface
  
      interface
        subroutine vlistDefVarZlevel(vlistID,varID,zlevel) bind(c,name='vlistDefVarZlevel')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: zlevel
       end subroutine vlistDefVarZlevel
      end interface
  
      interface
        integer(c_int) function vlistInqVarZlevel(vlistID,varID) bind(c,name='vlistInqVarZlevel')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarZlevel
      end interface
  
      interface
        subroutine vlistDefVarParam(vlistID,varID,param) bind(c,name='vlistDefVarParam')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: param
       end subroutine vlistDefVarParam
      end interface
  
      interface
        integer(c_int) function vlistInqVarParam(vlistID,varID) bind(c,name='vlistInqVarParam')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarParam
      end interface
  
      interface
        subroutine vlistDefVarCode(vlistID,varID,code) bind(c,name='vlistDefVarCode')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: code
       end subroutine vlistDefVarCode
      end interface
  
      interface
        integer(c_int) function vlistInqVarCode(vlistID,varID) bind(c,name='vlistInqVarCode')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarCode
      end interface
  
      interface
        subroutine vlistDefVarDatatype(vlistID,varID,datatype) bind(c,name='vlistDefVarDatatype')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: datatype
       end subroutine vlistDefVarDatatype
      end interface
  
      interface
        integer(c_int) function vlistInqVarDatatype(vlistID,varID) bind(c,name='vlistInqVarDatatype')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarDatatype
      end interface
  
      interface
        integer(c_int) function vlistInqVarNumber(vlistID,varID) bind(c,name='vlistInqVarNumber')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarNumber
      end interface
  
      interface
        subroutine vlistDefVarInstitut(vlistID,varID,instID) bind(c,name='vlistDefVarInstitut')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: instID
       end subroutine vlistDefVarInstitut
      end interface
  
      interface
        integer(c_int) function vlistInqVarInstitut(vlistID,varID) bind(c,name='vlistInqVarInstitut')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarInstitut
      end interface
  
      interface
        subroutine vlistDefVarModel(vlistID,varID,modelID) bind(c,name='vlistDefVarModel')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: modelID
       end subroutine vlistDefVarModel
      end interface
  
      interface
        integer(c_int) function vlistInqVarModel(vlistID,varID) bind(c,name='vlistInqVarModel')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarModel
      end interface
  
      interface
        subroutine vlistDefVarTable(vlistID,varID,tableID) bind(c,name='vlistDefVarTable')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: tableID
       end subroutine vlistDefVarTable
      end interface
  
      interface
        integer(c_int) function vlistInqVarTable(vlistID,varID) bind(c,name='vlistInqVarTable')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarTable
      end interface
  
      interface
        subroutine vlistDefVarName(vlistID,varID,name) bind(c,name='vlistDefVarName')
          import :: c_int,c_char
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: name
       end subroutine vlistDefVarName
      end interface
  
      interface
        subroutine vlistInqVarName(vlistID,varID,name) bind(c,name='vlistInqVarName')
          import :: c_int,c_char
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: name
       end subroutine vlistInqVarName
      end interface
  
      interface
        subroutine vlistDefVarStdname(vlistID,varID,stdname) bind(c,name='vlistDefVarStdname')
          import :: c_int,c_char
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: stdname
       end subroutine vlistDefVarStdname
      end interface
  
      interface
        subroutine vlistInqVarStdname(vlistID,varID,stdname) bind(c,name='vlistInqVarStdname')
          import :: c_int,c_char
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: stdname
       end subroutine vlistInqVarStdname
      end interface
  
      interface
        subroutine vlistDefVarLongname(vlistID,varID,longname) bind(c,name='vlistDefVarLongname')
          import :: c_int,c_char
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: longname
       end subroutine vlistDefVarLongname
      end interface
  
      interface
        subroutine vlistInqVarLongname(vlistID,varID,longname) bind(c,name='vlistInqVarLongname')
          import :: c_int,c_char
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: longname
       end subroutine vlistInqVarLongname
      end interface
  
      interface
        subroutine vlistDefVarUnits(vlistID,varID,units) bind(c,name='vlistDefVarUnits')
          import :: c_int,c_char
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: units
       end subroutine vlistDefVarUnits
      end interface
  
      interface
        subroutine vlistInqVarUnits(vlistID,varID,units) bind(c,name='vlistInqVarUnits')
          import :: c_int,c_char
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: units
       end subroutine vlistInqVarUnits
      end interface
  
      interface
        subroutine vlistDefVarMissval(vlistID,varID,missval) bind(c,name='vlistDefVarMissval')
          import :: c_int,c_double
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          real(c_double), value :: missval
       end subroutine vlistDefVarMissval
      end interface
  
      interface
        real(c_double) function vlistInqVarMissval(vlistID,varID) bind(c,name='vlistInqVarMissval')
          import :: c_int,c_double
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarMissval
      end interface
  
      interface
        subroutine vlistDefVarScalefactor(vlistID,varID,scalefactor) bind(c,name='vlistDefVarScalefactor')
          import :: c_int,c_double
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          real(c_double), value :: scalefactor
       end subroutine vlistDefVarScalefactor
      end interface
  
      interface
        real(c_double) function vlistInqVarScalefactor(vlistID,varID) bind(c,name='vlistInqVarScalefactor')
          import :: c_int,c_double
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarScalefactor
      end interface
  
      interface
        subroutine vlistDefVarAddoffset(vlistID,varID,addoffset) bind(c,name='vlistDefVarAddoffset')
          import :: c_int,c_double
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          real(c_double), value :: addoffset
       end subroutine vlistDefVarAddoffset
      end interface
  
      interface
        real(c_double) function vlistInqVarAddoffset(vlistID,varID) bind(c,name='vlistInqVarAddoffset')
          import :: c_int,c_double
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarAddoffset
      end interface
  
      interface
        subroutine vlistDefVarTsteptype(vlistID,varID,tsteptype) bind(c,name='vlistDefVarTsteptype')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: tsteptype
       end subroutine vlistDefVarTsteptype
      end interface
  
      interface
        integer(c_int) function vlistInqVarTsteptype(vlistID,varID) bind(c,name='vlistInqVarTsteptype')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarTsteptype
      end interface
  
      interface
        subroutine vlistDefVarTimave(vlistID,varID,timave) bind(c,name='vlistDefVarTimave')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: timave
       end subroutine vlistDefVarTimave
      end interface
  
      interface
        integer(c_int) function vlistInqVarTimave(vlistID,varID) bind(c,name='vlistInqVarTimave')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarTimave
      end interface
  
      interface
        subroutine vlistDefVarTimaccu(vlistID,varID,timaccu) bind(c,name='vlistDefVarTimaccu')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: timaccu
       end subroutine vlistDefVarTimaccu
      end interface
  
      interface
        integer(c_int) function vlistInqVarTimaccu(vlistID,varID) bind(c,name='vlistInqVarTimaccu')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarTimaccu
      end interface
  
      interface
        integer(c_int) function vlistInqVarSize(vlistID,varID) bind(c,name='vlistInqVarSize')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistInqVarSize
      end interface
  
      interface
        integer(c_int) function vlistInqVarID(vlistID,code) bind(c,name='vlistInqVarID')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: code
       end function vlistInqVarID
      end interface
  
      interface
        subroutine vlistDefIndex(vlistID,varID,levID,index) bind(c,name='vlistDefIndex')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: levID
          integer(c_int), value :: index
       end subroutine vlistDefIndex
      end interface
  
      interface
        integer(c_int) function vlistInqIndex(vlistID,varID,levID) bind(c,name='vlistInqIndex')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: levID
       end function vlistInqIndex
      end interface
  
      interface
        subroutine vlistDefFlag(vlistID,varID,levID,flag) bind(c,name='vlistDefFlag')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: levID
          integer(c_int), value :: flag
       end subroutine vlistDefFlag
      end interface
  
      interface
        integer(c_int) function vlistInqFlag(vlistID,varID,levID) bind(c,name='vlistInqFlag')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: levID
       end function vlistInqFlag
      end interface
  
      interface
        integer(c_int) function vlistFindVar(vlistID,fvarID) bind(c,name='vlistFindVar')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: fvarID
       end function vlistFindVar
      end interface
  
      interface
        integer(c_int) function vlistFindLevel(vlistID,fvarID,flevelID) bind(c,name='vlistFindLevel')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: fvarID
          integer(c_int), value :: flevelID
       end function vlistFindLevel
      end interface
  
      interface
        integer(c_int) function vlistMergedVar(vlistID,varID) bind(c,name='vlistMergedVar')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
       end function vlistMergedVar
      end interface
  
      interface
        integer(c_int) function vlistMergedLevel(vlistID,varID,levelID) bind(c,name='vlistMergedLevel')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: levelID
       end function vlistMergedLevel
      end interface
  
      interface
        integer(c_int) function vlistInqNatts(vlistID,varID,nattsp) bind(c,name='vlistInqNatts')
          import :: c_int
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), intent(out) :: nattsp
       end function vlistInqNatts
      end interface
  
      interface
        integer(c_int) function vlistInqAtt(vlistID,varID,attrnum,name,typep,lenp) bind(c,name='vlistInqAtt')
          import :: c_int,c_char
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          integer(c_int), value :: attrnum
          character(c_char), dimension(*) :: name
          integer(c_int), intent(out) :: typep
          integer(c_int), intent(out) :: lenp
       end function vlistInqAtt
      end interface
  
      interface
        integer(c_int) function vlistDelAtt(vlistID,varID,name) bind(c,name='vlistDelAtt')
          import :: c_int,c_char
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: name
       end function vlistDelAtt
      end interface
  
      interface
        integer(c_int) function vlistDefAttInt(vlistID,varID,name,len,ip_vec) bind(c,name='vlistDefAttInt')
          import :: c_int,c_char
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: name
          integer(c_int), value :: len
          integer(c_int), intent(in),dimension(*) :: ip_vec
       end function vlistDefAttInt
      end interface
  
      interface
        integer(c_int) function vlistDefAttFlt(vlistID,varID,name,len,dp_vec) bind(c,name='vlistDefAttFlt')
          import :: c_int,c_char,c_double
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: name
          integer(c_int), value :: len
          real(c_double), intent(in),dimension(*) :: dp_vec
       end function vlistDefAttFlt
      end interface
  
      interface
        integer(c_int) function vlistDefAttTxt(vlistID,varID,name,len,tp) bind(c,name='vlistDefAttTxt')
          import :: c_int,c_char
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: name
          integer(c_int), value :: len
          character(c_char), dimension(*) :: tp
       end function vlistDefAttTxt
      end interface
  
      interface
        integer(c_int) function vlistInqAttInt(vlistID,varID,name,mlen,ip_vec) bind(c,name='vlistInqAttInt')
          import :: c_int,c_char
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: name
          integer(c_int), value :: mlen
          integer(c_int), intent(out),dimension(*) :: ip_vec
       end function vlistInqAttInt
      end interface
  
      interface
        integer(c_int) function vlistInqAttFlt(vlistID,varID,name,mlen,dp_vec) bind(c,name='vlistInqAttFlt')
          import :: c_int,c_char,c_double
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: name
          integer(c_int), value :: mlen
          real(c_double), intent(out),dimension(*) :: dp_vec
       end function vlistInqAttFlt
      end interface
  
      interface
        integer(c_int) function vlistInqAttTxt(vlistID,varID,name,mlen,tp) bind(c,name='vlistInqAttTxt')
          import :: c_int,c_char
          integer(c_int), value :: vlistID
          integer(c_int), value :: varID
          character(c_char), dimension(*) :: name
          integer(c_int), value :: mlen
          character(c_char), dimension(*) :: tp
       end function vlistInqAttTxt
      end interface
  
      interface
        subroutine gridName(gridtype,gridnamev) bind(c,name='gridName')
          import :: c_int,c_char
          integer(c_int), value :: gridtype
          character(c_char), dimension(*) :: gridnamev
       end subroutine gridName
      end interface
  
      interface
        character(c_char) function gridNamePtr(gridtype) bind(c,name='gridNamePtr')
          import :: c_int,c_char
          integer(c_int), value :: gridtype
       end function gridNamePtr
      end interface
  
      interface
        subroutine gridCompress(gridID) bind(c,name='gridCompress')
          import :: c_int
          integer(c_int), value :: gridID
       end subroutine gridCompress
      end interface
  
      interface
        subroutine gridDefMask(gridID,mask_vec) bind(c,name='gridDefMask')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), intent(in),dimension(*) :: mask_vec
       end subroutine gridDefMask
      end interface
  
      interface
        integer(c_int) function gridInqMask(gridID,mask_vec) bind(c,name='gridInqMask')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), intent(out),dimension(*) :: mask_vec
       end function gridInqMask
      end interface
  
      interface
        subroutine gridPrint(gridID,opt) bind(c,name='gridPrint')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), value :: opt
       end subroutine gridPrint
      end interface
  
      interface
        integer(c_int) function gridSize() bind(c,name='gridSize')
          import :: c_int
       end function gridSize
      end interface
  
      interface
        integer(c_int) function gridCreate(gridtype,size) bind(c,name='gridCreate')
          import :: c_int
          integer(c_int), value :: gridtype
          integer(c_int), value :: size
       end function gridCreate
      end interface
  
      interface
        subroutine gridDestroy(gridID) bind(c,name='gridDestroy')
          import :: c_int
          integer(c_int), value :: gridID
       end subroutine gridDestroy
      end interface
  
      interface
        integer(c_int) function gridDuplicate(gridID) bind(c,name='gridDuplicate')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridDuplicate
      end interface
  
      interface
        integer(c_int) function gridInqType(gridID) bind(c,name='gridInqType')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridInqType
      end interface
  
      interface
        integer(c_int) function gridInqSize(gridID) bind(c,name='gridInqSize')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridInqSize
      end interface
  
      interface
        subroutine gridDefXsize(gridID,xsize) bind(c,name='gridDefXsize')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), value :: xsize
       end subroutine gridDefXsize
      end interface
  
      interface
        integer(c_int) function gridInqXsize(gridID) bind(c,name='gridInqXsize')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridInqXsize
      end interface
  
      interface
        subroutine gridDefYsize(gridID,ysize) bind(c,name='gridDefYsize')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), value :: ysize
       end subroutine gridDefYsize
      end interface
  
      interface
        integer(c_int) function gridInqYsize(gridID) bind(c,name='gridInqYsize')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridInqYsize
      end interface
  
      interface
        subroutine gridDefXvals(gridID,xvals_vec) bind(c,name='gridDefXvals')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), intent(in),dimension(*) :: xvals_vec
       end subroutine gridDefXvals
      end interface
  
      interface
        integer(c_int) function gridInqXvals(gridID,xvals_vec) bind(c,name='gridInqXvals')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), intent(out),dimension(*) :: xvals_vec
       end function gridInqXvals
      end interface
  
      interface
        subroutine gridDefYvals(gridID,yvals_vec) bind(c,name='gridDefYvals')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), intent(in),dimension(*) :: yvals_vec
       end subroutine gridDefYvals
      end interface
  
      interface
        integer(c_int) function gridInqYvals(gridID,yvals_vec) bind(c,name='gridInqYvals')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), intent(out),dimension(*) :: yvals_vec
       end function gridInqYvals
      end interface
  
      interface
        subroutine gridDefXname(gridID,xname) bind(c,name='gridDefXname')
          import :: c_int,c_char
          integer(c_int), value :: gridID
          character(c_char), dimension(*) :: xname
       end subroutine gridDefXname
      end interface
  
      interface
        subroutine gridDefXlongname(gridID,xlongname) bind(c,name='gridDefXlongname')
          import :: c_int,c_char
          integer(c_int), value :: gridID
          character(c_char), dimension(*) :: xlongname
       end subroutine gridDefXlongname
      end interface
  
      interface
        subroutine gridDefXunits(gridID,xunits) bind(c,name='gridDefXunits')
          import :: c_int,c_char
          integer(c_int), value :: gridID
          character(c_char), dimension(*) :: xunits
       end subroutine gridDefXunits
      end interface
  
      interface
        subroutine gridDefYname(gridID,yname) bind(c,name='gridDefYname')
          import :: c_int,c_char
          integer(c_int), value :: gridID
          character(c_char), dimension(*) :: yname
       end subroutine gridDefYname
      end interface
  
      interface
        subroutine gridDefYlongname(gridID,ylongname) bind(c,name='gridDefYlongname')
          import :: c_int,c_char
          integer(c_int), value :: gridID
          character(c_char), dimension(*) :: ylongname
       end subroutine gridDefYlongname
      end interface
  
      interface
        subroutine gridDefYunits(gridID,yunits) bind(c,name='gridDefYunits')
          import :: c_int,c_char
          integer(c_int), value :: gridID
          character(c_char), dimension(*) :: yunits
       end subroutine gridDefYunits
      end interface
  
      interface
        subroutine gridInqXname(gridID,xname) bind(c,name='gridInqXname')
          import :: c_int,c_char
          integer(c_int), value :: gridID
          character(c_char), dimension(*) :: xname
       end subroutine gridInqXname
      end interface
  
      interface
        subroutine gridInqXlongname(gridID,xlongname) bind(c,name='gridInqXlongname')
          import :: c_int,c_char
          integer(c_int), value :: gridID
          character(c_char), dimension(*) :: xlongname
       end subroutine gridInqXlongname
      end interface
  
      interface
        subroutine gridInqXstdname(gridID,xstdname) bind(c,name='gridInqXstdname')
          import :: c_int,c_char
          integer(c_int), value :: gridID
          character(c_char), dimension(*) :: xstdname
       end subroutine gridInqXstdname
      end interface
  
      interface
        subroutine gridInqXunits(gridID,xunits) bind(c,name='gridInqXunits')
          import :: c_int,c_char
          integer(c_int), value :: gridID
          character(c_char), dimension(*) :: xunits
       end subroutine gridInqXunits
      end interface
  
      interface
        subroutine gridInqYname(gridID,yname) bind(c,name='gridInqYname')
          import :: c_int,c_char
          integer(c_int), value :: gridID
          character(c_char), dimension(*) :: yname
       end subroutine gridInqYname
      end interface
  
      interface
        subroutine gridInqYlongname(gridID,ylongname) bind(c,name='gridInqYlongname')
          import :: c_int,c_char
          integer(c_int), value :: gridID
          character(c_char), dimension(*) :: ylongname
       end subroutine gridInqYlongname
      end interface
  
      interface
        subroutine gridInqYstdname(gridID,ystdname) bind(c,name='gridInqYstdname')
          import :: c_int,c_char
          integer(c_int), value :: gridID
          character(c_char), dimension(*) :: ystdname
       end subroutine gridInqYstdname
      end interface
  
      interface
        subroutine gridInqYunits(gridID,yunits) bind(c,name='gridInqYunits')
          import :: c_int,c_char
          integer(c_int), value :: gridID
          character(c_char), dimension(*) :: yunits
       end subroutine gridInqYunits
      end interface
  
      interface
        subroutine gridDefPrec(gridID,prec) bind(c,name='gridDefPrec')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), value :: prec
       end subroutine gridDefPrec
      end interface
  
      interface
        integer(c_int) function gridInqPrec(gridID) bind(c,name='gridInqPrec')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridInqPrec
      end interface
  
      interface
        real(c_double) function gridInqXval(gridID,index) bind(c,name='gridInqXval')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          integer(c_int), value :: index
       end function gridInqXval
      end interface
  
      interface
        real(c_double) function gridInqYval(gridID,index) bind(c,name='gridInqYval')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          integer(c_int), value :: index
       end function gridInqYval
      end interface
  
      interface
        real(c_double) function gridInqXinc(gridID) bind(c,name='gridInqXinc')
          import :: c_int,c_double
          integer(c_int), value :: gridID
       end function gridInqXinc
      end interface
  
      interface
        real(c_double) function gridInqYinc(gridID) bind(c,name='gridInqYinc')
          import :: c_int,c_double
          integer(c_int), value :: gridID
       end function gridInqYinc
      end interface
  
      interface
        integer(c_int) function gridIsCircular(gridID) bind(c,name='gridIsCircular')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridIsCircular
      end interface
  
      interface
        integer(c_int) function gridIsRotated(gridID) bind(c,name='gridIsRotated')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridIsRotated
      end interface
  
      interface
        real(c_double) function gridInqXpole(gridID) bind(c,name='gridInqXpole')
          import :: c_int,c_double
          integer(c_int), value :: gridID
       end function gridInqXpole
      end interface
  
      interface
        subroutine gridDefXpole(gridID,xpole) bind(c,name='gridDefXpole')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), value :: xpole
       end subroutine gridDefXpole
      end interface
  
      interface
        real(c_double) function gridInqYpole(gridID) bind(c,name='gridInqYpole')
          import :: c_int,c_double
          integer(c_int), value :: gridID
       end function gridInqYpole
      end interface
  
      interface
        subroutine gridDefYpole(gridID,ypole) bind(c,name='gridDefYpole')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), value :: ypole
       end subroutine gridDefYpole
      end interface
  
      interface
        real(c_double) function gridInqAngle(gridID) bind(c,name='gridInqAngle')
          import :: c_int,c_double
          integer(c_int), value :: gridID
       end function gridInqAngle
      end interface
  
      interface
        subroutine gridDefAngle(gridID,angle) bind(c,name='gridDefAngle')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), value :: angle
       end subroutine gridDefAngle
      end interface
  
      interface
        subroutine gridDefTrunc(gridID,trunc) bind(c,name='gridDefTrunc')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), value :: trunc
       end subroutine gridDefTrunc
      end interface
  
      interface
        integer(c_int) function gridInqTrunc(gridID) bind(c,name='gridInqTrunc')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridInqTrunc
      end interface
  
      interface
        integer(c_int) function gridInqGMEnd(gridID) bind(c,name='gridInqGMEnd')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridInqGMEnd
      end interface
  
      interface
        subroutine gridDefGMEnd(gridID,nd) bind(c,name='gridDefGMEnd')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), value :: nd
       end subroutine gridDefGMEnd
      end interface
  
      interface
        integer(c_int) function gridInqGMEni(gridID) bind(c,name='gridInqGMEni')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridInqGMEni
      end interface
  
      interface
        subroutine gridDefGMEni(gridID,ni) bind(c,name='gridDefGMEni')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), value :: ni
       end subroutine gridDefGMEni
      end interface
  
      interface
        integer(c_int) function gridInqGMEni2(gridID) bind(c,name='gridInqGMEni2')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridInqGMEni2
      end interface
  
      interface
        subroutine gridDefGMEni2(gridID,ni2) bind(c,name='gridDefGMEni2')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), value :: ni2
       end subroutine gridDefGMEni2
      end interface
  
      interface
        integer(c_int) function gridInqGMEni3(gridID) bind(c,name='gridInqGMEni3')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridInqGMEni3
      end interface
  
      interface
        subroutine gridDefGMEni3(gridID,ni3) bind(c,name='gridDefGMEni3')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), value :: ni3
       end subroutine gridDefGMEni3
      end interface
  
      interface
        subroutine gridDefLCC(gridID,originLon,originLat,lonParY,lat1,lat2,xinc,yinc,projflag,scanflag) bind(c,name='gridDefLCC')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), value :: originLon
          real(c_double), value :: originLat
          real(c_double), value :: lonParY
          real(c_double), value :: lat1
          real(c_double), value :: lat2
          real(c_double), value :: xinc
          real(c_double), value :: yinc
          integer(c_int), value :: projflag
          integer(c_int), value :: scanflag
       end subroutine gridDefLCC
      end interface
  
      interface
        subroutine gridInqLCC(gridID,originLon,originLat,lonParY,lat1,lat2,xinc,yinc,projflag,scanflag) bind(c,name='gridInqLCC')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), intent(out) :: originLon
          real(c_double), intent(out) :: originLat
          real(c_double), intent(out) :: lonParY
          real(c_double), intent(out) :: lat1
          real(c_double), intent(out) :: lat2
          real(c_double), intent(out) :: xinc
          real(c_double), intent(out) :: yinc
          integer(c_int), intent(out) :: projflag
          integer(c_int), intent(out) :: scanflag
       end subroutine gridInqLCC
      end interface
  
      interface
        subroutine gridDefLcc2(gridID,earth_radius,lon_0,lat_0,lat_1,lat_2) bind(c,name='gridDefLcc2')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), value :: earth_radius
          real(c_double), value :: lon_0
          real(c_double), value :: lat_0
          real(c_double), value :: lat_1
          real(c_double), value :: lat_2
       end subroutine gridDefLcc2
      end interface
  
      interface
        subroutine gridInqLcc2(gridID,earth_radius,lon_0,lat_0,lat_1,lat_2) bind(c,name='gridInqLcc2')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), intent(out) :: earth_radius
          real(c_double), intent(out) :: lon_0
          real(c_double), intent(out) :: lat_0
          real(c_double), intent(out) :: lat_1
          real(c_double), intent(out) :: lat_2
       end subroutine gridInqLcc2
      end interface
  
      interface
        subroutine gridDefLaea(gridID,earth_radius,lon_0,lat_0) bind(c,name='gridDefLaea')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), value :: earth_radius
          real(c_double), value :: lon_0
          real(c_double), value :: lat_0
       end subroutine gridDefLaea
      end interface
  
      interface
        subroutine gridInqLaea(gridID,earth_radius,lon_0,lat_0) bind(c,name='gridInqLaea')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), intent(out) :: earth_radius
          real(c_double), intent(out) :: lon_0
          real(c_double), intent(out) :: lat_0
       end subroutine gridInqLaea
      end interface
  
      interface
        subroutine gridDefArea(gridID,area_vec) bind(c,name='gridDefArea')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), intent(in),dimension(*) :: area_vec
       end subroutine gridDefArea
      end interface
  
      interface
        subroutine gridInqArea(gridID,area_vec) bind(c,name='gridInqArea')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), intent(out),dimension(*) :: area_vec
       end subroutine gridInqArea
      end interface
  
      interface
        integer(c_int) function gridHasArea(gridID) bind(c,name='gridHasArea')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridHasArea
      end interface
  
      interface
        subroutine gridDefNvertex(gridID,nvertex) bind(c,name='gridDefNvertex')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), value :: nvertex
       end subroutine gridDefNvertex
      end interface
  
      interface
        integer(c_int) function gridInqNvertex(gridID) bind(c,name='gridInqNvertex')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridInqNvertex
      end interface
  
      interface
        subroutine gridDefXbounds(gridID,xbounds_vec) bind(c,name='gridDefXbounds')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), intent(in),dimension(*) :: xbounds_vec
       end subroutine gridDefXbounds
      end interface
  
      interface
        integer(c_int) function gridInqXbounds(gridID,xbounds_vec) bind(c,name='gridInqXbounds')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), intent(out),dimension(*) :: xbounds_vec
       end function gridInqXbounds
      end interface
  
      interface
        subroutine gridDefYbounds(gridID,ybounds_vec) bind(c,name='gridDefYbounds')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), intent(in),dimension(*) :: ybounds_vec
       end subroutine gridDefYbounds
      end interface
  
      interface
        integer(c_int) function gridInqYbounds(gridID,ybounds_vec) bind(c,name='gridInqYbounds')
          import :: c_int,c_double
          integer(c_int), value :: gridID
          real(c_double), intent(out),dimension(*) :: ybounds_vec
       end function gridInqYbounds
      end interface
  
      interface
        subroutine gridDefRowlon(gridID,nrowlon,rowlon_vec) bind(c,name='gridDefRowlon')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), value :: nrowlon
          integer(c_int), intent(in),dimension(*) :: rowlon_vec
       end subroutine gridDefRowlon
      end interface
  
      interface
        subroutine gridInqRowlon(gridID,rowlon_vec) bind(c,name='gridInqRowlon')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), intent(out),dimension(*) :: rowlon_vec
       end subroutine gridInqRowlon
      end interface
  
      interface
        subroutine gridChangeType(gridID,gridtype) bind(c,name='gridChangeType')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), value :: gridtype
       end subroutine gridChangeType
      end interface
  
      interface
        subroutine gridDefComplexPacking(gridID,lpack) bind(c,name='gridDefComplexPacking')
          import :: c_int
          integer(c_int), value :: gridID
          integer(c_int), value :: lpack
       end subroutine gridDefComplexPacking
      end interface
  
      interface
        integer(c_int) function gridInqComplexPacking(gridID) bind(c,name='gridInqComplexPacking')
          import :: c_int
          integer(c_int), value :: gridID
       end function gridInqComplexPacking
      end interface
  
      interface
        subroutine zaxisName(zaxistype,zaxisnamev) bind(c,name='zaxisName')
          import :: c_int,c_char
          integer(c_int), value :: zaxistype
          character(c_char), dimension(*) :: zaxisnamev
       end subroutine zaxisName
      end interface
  
      interface
        integer(c_int) function zaxisCreate(zaxistype,size) bind(c,name='zaxisCreate')
          import :: c_int
          integer(c_int), value :: zaxistype
          integer(c_int), value :: size
       end function zaxisCreate
      end interface
  
      interface
        subroutine zaxisDestroy(zaxisID) bind(c,name='zaxisDestroy')
          import :: c_int
          integer(c_int), value :: zaxisID
       end subroutine zaxisDestroy
      end interface
  
      interface
        integer(c_int) function zaxisInqType(zaxisID) bind(c,name='zaxisInqType')
          import :: c_int
          integer(c_int), value :: zaxisID
       end function zaxisInqType
      end interface
  
      interface
        integer(c_int) function zaxisInqSize(zaxisID) bind(c,name='zaxisInqSize')
          import :: c_int
          integer(c_int), value :: zaxisID
       end function zaxisInqSize
      end interface
  
      interface
        integer(c_int) function zaxisDuplicate(zaxisID) bind(c,name='zaxisDuplicate')
          import :: c_int
          integer(c_int), value :: zaxisID
       end function zaxisDuplicate
      end interface
  
      interface
        subroutine zaxisResize(zaxisID,size) bind(c,name='zaxisResize')
          import :: c_int
          integer(c_int), value :: zaxisID
          integer(c_int), value :: size
       end subroutine zaxisResize
      end interface
  
      interface
        subroutine zaxisPrint(zaxisID) bind(c,name='zaxisPrint')
          import :: c_int
          integer(c_int), value :: zaxisID
       end subroutine zaxisPrint
      end interface
  
      interface
        integer(c_int) function zaxisSize() bind(c,name='zaxisSize')
          import :: c_int
       end function zaxisSize
      end interface
  
      interface
        subroutine zaxisDefLevels(zaxisID,levels_vec) bind(c,name='zaxisDefLevels')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
          real(c_double), intent(in),dimension(*) :: levels_vec
       end subroutine zaxisDefLevels
      end interface
  
      interface
        subroutine zaxisInqLevels(zaxisID,levels_vec) bind(c,name='zaxisInqLevels')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
          real(c_double), intent(out),dimension(*) :: levels_vec
       end subroutine zaxisInqLevels
      end interface
  
      interface
        subroutine zaxisDefLevel(zaxisID,levelID,levels) bind(c,name='zaxisDefLevel')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
          integer(c_int), value :: levelID
          real(c_double), value :: levels
       end subroutine zaxisDefLevel
      end interface
  
      interface
        real(c_double) function zaxisInqLevel(zaxisID,levelID) bind(c,name='zaxisInqLevel')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
          integer(c_int), value :: levelID
       end function zaxisInqLevel
      end interface
  
      interface
        subroutine zaxisDefName(zaxisID,name) bind(c,name='zaxisDefName')
          import :: c_int,c_char
          integer(c_int), value :: zaxisID
          character(c_char), dimension(*) :: name
       end subroutine zaxisDefName
      end interface
  
      interface
        subroutine zaxisDefLongname(zaxisID,longname) bind(c,name='zaxisDefLongname')
          import :: c_int,c_char
          integer(c_int), value :: zaxisID
          character(c_char), dimension(*) :: longname
       end subroutine zaxisDefLongname
      end interface
  
      interface
        subroutine zaxisDefUnits(zaxisID,units) bind(c,name='zaxisDefUnits')
          import :: c_int,c_char
          integer(c_int), value :: zaxisID
          character(c_char), dimension(*) :: units
       end subroutine zaxisDefUnits
      end interface
  
      interface
        subroutine zaxisInqName(zaxisID,name) bind(c,name='zaxisInqName')
          import :: c_int,c_char
          integer(c_int), value :: zaxisID
          character(c_char), dimension(*) :: name
       end subroutine zaxisInqName
      end interface
  
      interface
        subroutine zaxisInqLongname(zaxisID,longname) bind(c,name='zaxisInqLongname')
          import :: c_int,c_char
          integer(c_int), value :: zaxisID
          character(c_char), dimension(*) :: longname
       end subroutine zaxisInqLongname
      end interface
  
      interface
        subroutine zaxisInqUnits(zaxisID,units) bind(c,name='zaxisInqUnits')
          import :: c_int,c_char
          integer(c_int), value :: zaxisID
          character(c_char), dimension(*) :: units
       end subroutine zaxisInqUnits
      end interface
  
      interface
        subroutine zaxisDefPrec(zaxisID,prec) bind(c,name='zaxisDefPrec')
          import :: c_int
          integer(c_int), value :: zaxisID
          integer(c_int), value :: prec
       end subroutine zaxisDefPrec
      end interface
  
      interface
        integer(c_int) function zaxisInqPrec(zaxisID) bind(c,name='zaxisInqPrec')
          import :: c_int
          integer(c_int), value :: zaxisID
       end function zaxisInqPrec
      end interface
  
      interface
        subroutine zaxisDefLtype(zaxisID,ltype) bind(c,name='zaxisDefLtype')
          import :: c_int
          integer(c_int), value :: zaxisID
          integer(c_int), value :: ltype
       end subroutine zaxisDefLtype
      end interface
  
      interface
        integer(c_int) function zaxisInqLtype(zaxisID) bind(c,name='zaxisInqLtype')
          import :: c_int
          integer(c_int), value :: zaxisID
       end function zaxisInqLtype
      end interface
  
      interface
        real(c_double) function zaxisInqLevelsPtr(zaxisID) bind(c,name='zaxisInqLevelsPtr')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
       end function zaxisInqLevelsPtr
      end interface
  
      interface
        subroutine zaxisDefVct(zaxisID,size,vct_vec) bind(c,name='zaxisDefVct')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
          integer(c_int), value :: size
          real(c_double), intent(in),dimension(*) :: vct_vec
       end subroutine zaxisDefVct
      end interface
  
      interface
        integer(c_int) function zaxisInqVctSize(zaxisID) bind(c,name='zaxisInqVctSize')
          import :: c_int
          integer(c_int), value :: zaxisID
       end function zaxisInqVctSize
      end interface
  
      interface
        real(c_double) function zaxisInqVctPtr(zaxisID) bind(c,name='zaxisInqVctPtr')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
       end function zaxisInqVctPtr
      end interface
  
      interface
        integer(c_int) function zaxisInqLbounds(zaxisID,lbounds_vec) bind(c,name='zaxisInqLbounds')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
          real(c_double), intent(out),dimension(*) :: lbounds_vec
       end function zaxisInqLbounds
      end interface
  
      interface
        integer(c_int) function zaxisInqUbounds(zaxisID,ubounds_vec) bind(c,name='zaxisInqUbounds')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
          real(c_double), intent(out),dimension(*) :: ubounds_vec
       end function zaxisInqUbounds
      end interface
  
      interface
        integer(c_int) function zaxisInqWeights(zaxisID,weights_vec) bind(c,name='zaxisInqWeights')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
          real(c_double), intent(out),dimension(*) :: weights_vec
       end function zaxisInqWeights
      end interface
  
      interface
        real(c_double) function zaxisInqLbound(zaxisID,index) bind(c,name='zaxisInqLbound')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
          integer(c_int), value :: index
       end function zaxisInqLbound
      end interface
  
      interface
        real(c_double) function zaxisInqUbound(zaxisID,index) bind(c,name='zaxisInqUbound')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
          integer(c_int), value :: index
       end function zaxisInqUbound
      end interface
  
      interface
        subroutine zaxisDefLbounds(zaxisID,lbounds_vec) bind(c,name='zaxisDefLbounds')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
          real(c_double), intent(in),dimension(*) :: lbounds_vec
       end subroutine zaxisDefLbounds
      end interface
  
      interface
        subroutine zaxisDefUbounds(zaxisID,ubounds_vec) bind(c,name='zaxisDefUbounds')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
          real(c_double), intent(in),dimension(*) :: ubounds_vec
       end subroutine zaxisDefUbounds
      end interface
  
      interface
        subroutine zaxisDefWeights(zaxisID,weights_vec) bind(c,name='zaxisDefWeights')
          import :: c_int,c_double
          integer(c_int), value :: zaxisID
          real(c_double), intent(in),dimension(*) :: weights_vec
       end subroutine zaxisDefWeights
      end interface
  
      interface
        subroutine zaxisChangeType(zaxisID,zaxistype) bind(c,name='zaxisChangeType')
          import :: c_int
          integer(c_int), value :: zaxisID
          integer(c_int), value :: zaxistype
       end subroutine zaxisChangeType
      end interface
  
      interface
        integer(c_int) function taxisCreate(timetype) bind(c,name='taxisCreate')
          import :: c_int
          integer(c_int), value :: timetype
       end function taxisCreate
      end interface
  
      interface
        subroutine taxisDestroy(taxisID) bind(c,name='taxisDestroy')
          import :: c_int
          integer(c_int), value :: taxisID
       end subroutine taxisDestroy
      end interface
  
      interface
        integer(c_int) function taxisDuplicate(taxisID) bind(c,name='taxisDuplicate')
          import :: c_int
          integer(c_int), value :: taxisID
       end function taxisDuplicate
      end interface
  
      interface
        subroutine taxisCopyTimestep(taxisIDdes,taxisIDsrc) bind(c,name='taxisCopyTimestep')
          import :: c_int
          integer(c_int), value :: taxisIDdes
          integer(c_int), value :: taxisIDsrc
       end subroutine taxisCopyTimestep
      end interface
  
      interface
        subroutine taxisDefType(taxisID,type) bind(c,name='taxisDefType')
          import :: c_int
          integer(c_int), value :: taxisID
          integer(c_int), value :: type
       end subroutine taxisDefType
      end interface
  
      interface
        subroutine taxisDefVdate(taxisID,date) bind(c,name='taxisDefVdate')
          import :: c_int
          integer(c_int), value :: taxisID
          integer(c_int), value :: date
       end subroutine taxisDefVdate
      end interface
  
      interface
        subroutine taxisDefVtime(taxisID,time) bind(c,name='taxisDefVtime')
          import :: c_int
          integer(c_int), value :: taxisID
          integer(c_int), value :: time
       end subroutine taxisDefVtime
      end interface
  
      interface
        subroutine taxisDefRdate(taxisID,date) bind(c,name='taxisDefRdate')
          import :: c_int
          integer(c_int), value :: taxisID
          integer(c_int), value :: date
       end subroutine taxisDefRdate
      end interface
  
      interface
        subroutine taxisDefRtime(taxisID,time) bind(c,name='taxisDefRtime')
          import :: c_int
          integer(c_int), value :: taxisID
          integer(c_int), value :: time
       end subroutine taxisDefRtime
      end interface
  
      interface
        integer(c_int) function taxisHasBounds(taxisID) bind(c,name='taxisHasBounds')
          import :: c_int
          integer(c_int), value :: taxisID
       end function taxisHasBounds
      end interface
  
      interface
        subroutine taxisDeleteBounds(taxisID) bind(c,name='taxisDeleteBounds')
          import :: c_int
          integer(c_int), value :: taxisID
       end subroutine taxisDeleteBounds
      end interface
  
      interface
        subroutine taxisDefVdateBounds(taxisID,vdate_lb,vdate_ub) bind(c,name='taxisDefVdateBounds')
          import :: c_int
          integer(c_int), value :: taxisID
          integer(c_int), value :: vdate_lb
          integer(c_int), value :: vdate_ub
       end subroutine taxisDefVdateBounds
      end interface
  
      interface
        subroutine taxisDefVtimeBounds(taxisID,vtime_lb,vtime_ub) bind(c,name='taxisDefVtimeBounds')
          import :: c_int
          integer(c_int), value :: taxisID
          integer(c_int), value :: vtime_lb
          integer(c_int), value :: vtime_ub
       end subroutine taxisDefVtimeBounds
      end interface
  
      interface
        subroutine taxisInqVdateBounds(taxisID,vdate_lb,vdate_ub) bind(c,name='taxisInqVdateBounds')
          import :: c_int
          integer(c_int), value :: taxisID
          integer(c_int), intent(out) :: vdate_lb
          integer(c_int), intent(out) :: vdate_ub
       end subroutine taxisInqVdateBounds
      end interface
  
      interface
        subroutine taxisInqVtimeBounds(taxisID,vtime_lb,vtime_ub) bind(c,name='taxisInqVtimeBounds')
          import :: c_int
          integer(c_int), value :: taxisID
          integer(c_int), intent(out) :: vtime_lb
          integer(c_int), intent(out) :: vtime_ub
       end subroutine taxisInqVtimeBounds
      end interface
  
      interface
        subroutine taxisDefCalendar(taxisID,calendar) bind(c,name='taxisDefCalendar')
          import :: c_int
          integer(c_int), value :: taxisID
          integer(c_int), value :: calendar
       end subroutine taxisDefCalendar
      end interface
  
      interface
        subroutine taxisDefTunit(taxisID,tunit) bind(c,name='taxisDefTunit')
          import :: c_int
          integer(c_int), value :: taxisID
          integer(c_int), value :: tunit
       end subroutine taxisDefTunit
      end interface
  
      interface
        subroutine taxisDefNumavg(taxisID,numavg) bind(c,name='taxisDefNumavg')
          import :: c_int
          integer(c_int), value :: taxisID
          integer(c_int), value :: numavg
       end subroutine taxisDefNumavg
      end interface
  
      interface
        integer(c_int) function taxisInqType(taxisID) bind(c,name='taxisInqType')
          import :: c_int
          integer(c_int), value :: taxisID
       end function taxisInqType
      end interface
  
      interface
        integer(c_int) function taxisInqVdate(taxisID) bind(c,name='taxisInqVdate')
          import :: c_int
          integer(c_int), value :: taxisID
       end function taxisInqVdate
      end interface
  
      interface
        integer(c_int) function taxisInqVtime(taxisID) bind(c,name='taxisInqVtime')
          import :: c_int
          integer(c_int), value :: taxisID
       end function taxisInqVtime
      end interface
  
      interface
        integer(c_int) function taxisInqRdate(taxisID) bind(c,name='taxisInqRdate')
          import :: c_int
          integer(c_int), value :: taxisID
       end function taxisInqRdate
      end interface
  
      interface
        integer(c_int) function taxisInqRtime(taxisID) bind(c,name='taxisInqRtime')
          import :: c_int
          integer(c_int), value :: taxisID
       end function taxisInqRtime
      end interface
  
      interface
        integer(c_int) function taxisInqCalendar(taxisID) bind(c,name='taxisInqCalendar')
          import :: c_int
          integer(c_int), value :: taxisID
       end function taxisInqCalendar
      end interface
  
      interface
        integer(c_int) function taxisInqTunit(taxisID) bind(c,name='taxisInqTunit')
          import :: c_int
          integer(c_int), value :: taxisID
       end function taxisInqTunit
      end interface
  
      interface
        integer(c_int) function taxisInqNumavg(taxisID) bind(c,name='taxisInqNumavg')
          import :: c_int
          integer(c_int), value :: taxisID
       end function taxisInqNumavg
      end interface
  
      interface
        character(c_char) function tunitNamePtr(tunitID) bind(c,name='tunitNamePtr')
          import :: c_int,c_char
          integer(c_int), value :: tunitID
       end function tunitNamePtr
      end interface
  
      interface
        integer(c_int) function institutDef(center,subcenter,name,longname) bind(c,name='institutDef')
          import :: c_int,c_char
          integer(c_int), value :: center
          integer(c_int), value :: subcenter
          character(c_char), dimension(*) :: name
          character(c_char), dimension(*) :: longname
       end function institutDef
      end interface
  
      interface
        integer(c_int) function institutInq(center,subcenter,name,longname) bind(c,name='institutInq')
          import :: c_int,c_char
          integer(c_int), value :: center
          integer(c_int), value :: subcenter
          character(c_char), dimension(*) :: name
          character(c_char), dimension(*) :: longname
       end function institutInq
      end interface
  
      interface
        integer(c_int) function institutInqNumber() bind(c,name='institutInqNumber')
          import :: c_int
       end function institutInqNumber
      end interface
  
      interface
        integer(c_int) function institutInqCenter(instID) bind(c,name='institutInqCenter')
          import :: c_int
          integer(c_int), value :: instID
       end function institutInqCenter
      end interface
  
      interface
        integer(c_int) function institutInqSubcenter(instID) bind(c,name='institutInqSubcenter')
          import :: c_int
          integer(c_int), value :: instID
       end function institutInqSubcenter
      end interface
  
      interface
        character(c_char) function institutInqNamePtr(instID) bind(c,name='institutInqNamePtr')
          import :: c_int,c_char
          integer(c_int), value :: instID
       end function institutInqNamePtr
      end interface
  
      interface
        character(c_char) function institutInqLongnamePtr(instID) bind(c,name='institutInqLongnamePtr')
          import :: c_int,c_char
          integer(c_int), value :: instID
       end function institutInqLongnamePtr
      end interface
  
      interface
        integer(c_int) function modelDef(instID,modelgribID,name) bind(c,name='modelDef')
          import :: c_int,c_char
          integer(c_int), value :: instID
          integer(c_int), value :: modelgribID
          character(c_char), dimension(*) :: name
       end function modelDef
      end interface
  
      interface
        integer(c_int) function modelInq(instID,modelgribID,name) bind(c,name='modelInq')
          import :: c_int,c_char
          integer(c_int), value :: instID
          integer(c_int), value :: modelgribID
          character(c_char), dimension(*) :: name
       end function modelInq
      end interface
  
      interface
        integer(c_int) function modelInqInstitut(modelID) bind(c,name='modelInqInstitut')
          import :: c_int
          integer(c_int), value :: modelID
       end function modelInqInstitut
      end interface
  
      interface
        integer(c_int) function modelInqGribID(modelID) bind(c,name='modelInqGribID')
          import :: c_int
          integer(c_int), value :: modelID
       end function modelInqGribID
      end interface
  
      interface
        character(c_char) function modelInqNamePtr(modelID) bind(c,name='modelInqNamePtr')
          import :: c_int,c_char
          integer(c_int), value :: modelID
       end function modelInqNamePtr
      end interface
  
      interface
        subroutine tableWriteC(filename,tableID) bind(c,name='tableWriteC')
          import :: c_char,c_int
          character(c_char), dimension(*) :: filename
          integer(c_int), value :: tableID
       end subroutine tableWriteC
      end interface
  
      interface
        subroutine tableWrite(filename,tableID) bind(c,name='tableWrite')
          import :: c_char,c_int
          character(c_char), dimension(*) :: filename
          integer(c_int), value :: tableID
       end subroutine tableWrite
      end interface
  
      interface
        integer(c_int) function tableRead(tablefile) bind(c,name='tableRead')
          import :: c_char,c_int
          character(c_char), dimension(*) :: tablefile
       end function tableRead
      end interface
  
      interface
        integer(c_int) function tableDef(modelID,tablenum,tablename) bind(c,name='tableDef')
          import :: c_int,c_char
          integer(c_int), value :: modelID
          integer(c_int), value :: tablenum
          character(c_char), dimension(*) :: tablename
       end function tableDef
      end interface
  
      interface
        character(c_char) function tableInqNamePtr(tableID) bind(c,name='tableInqNamePtr')
          import :: c_int,c_char
          integer(c_int), value :: tableID
       end function tableInqNamePtr
      end interface
  
      interface
        subroutine tableDefEntry(tableID,code,name,longname,units) bind(c,name='tableDefEntry')
          import :: c_int,c_char
          integer(c_int), value :: tableID
          integer(c_int), value :: code
          character(c_char), dimension(*) :: name
          character(c_char), dimension(*) :: longname
          character(c_char), dimension(*) :: units
       end subroutine tableDefEntry
      end interface
  
      interface
        integer(c_int) function tableInq(modelID,tablenum,tablename) bind(c,name='tableInq')
          import :: c_int,c_char
          integer(c_int), value :: modelID
          integer(c_int), value :: tablenum
          character(c_char), dimension(*) :: tablename
       end function tableInq
      end interface
  
      interface
        integer(c_int) function tableInqNumber() bind(c,name='tableInqNumber')
          import :: c_int
       end function tableInqNumber
      end interface
  
      interface
        integer(c_int) function tableInqNum(tableID) bind(c,name='tableInqNum')
          import :: c_int
          integer(c_int), value :: tableID
       end function tableInqNum
      end interface
  
      interface
        integer(c_int) function tableInqModel(tableID) bind(c,name='tableInqModel')
          import :: c_int
          integer(c_int), value :: tableID
       end function tableInqModel
      end interface
  
      interface
        subroutine tableInqPar(tableID,code,name,longname,units) bind(c,name='tableInqPar')
          import :: c_int,c_char
          integer(c_int), value :: tableID
          integer(c_int), value :: code
          character(c_char), dimension(*) :: name
          character(c_char), dimension(*) :: longname
          character(c_char), dimension(*) :: units
       end subroutine tableInqPar
      end interface
  
      interface
        integer(c_int) function tableInqParCode(tableID,name,code) bind(c,name='tableInqParCode')
          import :: c_int,c_char
          integer(c_int), value :: tableID
          character(c_char), dimension(*) :: name
          integer(c_int), intent(out) :: code
       end function tableInqParCode
      end interface
  
      interface
        integer(c_int) function tableInqParName(tableID,code,name) bind(c,name='tableInqParName')
          import :: c_int,c_char
          integer(c_int), value :: tableID
          integer(c_int), value :: code
          character(c_char), dimension(*) :: name
       end function tableInqParName
      end interface
  
      interface
        integer(c_int) function tableInqParLongname(tableID,code,longname) bind(c,name='tableInqParLongname')
          import :: c_int,c_char
          integer(c_int), value :: tableID
          integer(c_int), value :: code
          character(c_char), dimension(*) :: longname
       end function tableInqParLongname
      end interface
  
      interface
        integer(c_int) function tableInqParUnits(tableID,code,units) bind(c,name='tableInqParUnits')
          import :: c_int,c_char
          integer(c_int), value :: tableID
          integer(c_int), value :: code
          character(c_char), dimension(*) :: units
       end function tableInqParUnits
      end interface
  
      interface
        character(c_char) function tableInqParNamePtr(tableID,parID) bind(c,name='tableInqParNamePtr')
          import :: c_int,c_char
          integer(c_int), value :: tableID
          integer(c_int), value :: parID
       end function tableInqParNamePtr
      end interface
  
      interface
        character(c_char) function tableInqParLongnamePtr(tableID,parID) bind(c,name='tableInqParLongnamePtr')
          import :: c_int,c_char
          integer(c_int), value :: tableID
          integer(c_int), value :: parID
       end function tableInqParLongnamePtr
      end interface
  
      interface
        character(c_char) function tableInqParUnitsPtr(tableID,parID) bind(c,name='tableInqParUnitsPtr')
          import :: c_int,c_char
          integer(c_int), value :: tableID
          integer(c_int), value :: parID
       end function tableInqParUnitsPtr
      end interface
  
      interface
        subroutine streamDefHistory(streamID,size,history) bind(c,name='streamDefHistory')
          import :: c_int,c_char
          integer(c_int), value :: streamID
          integer(c_int), value :: size
          character(c_char), dimension(*) :: history
       end subroutine streamDefHistory
      end interface
  
      interface
        integer(c_int) function streamInqHistorySize(streamID) bind(c,name='streamInqHistorySize')
          import :: c_int
          integer(c_int), value :: streamID
       end function streamInqHistorySize
      end interface
  
      interface
        subroutine streamInqHistoryString(streamID,history) bind(c,name='streamInqHistoryString')
          import :: c_int,c_char
          integer(c_int), value :: streamID
          character(c_char), dimension(*) :: history
       end subroutine streamInqHistoryString
      end interface
  
      public :: cdiStringError
      public :: cdiDebug
      public :: cdiLibraryVersion
      public :: cdiPrintVersion
      public :: cdiDefMissval
      public :: cdiInqMissval
      public :: cdiDefGlobal
      public :: cdiParamToString
      public :: cdiDecodeParam
      public :: cdiEncodeParam
      public :: cdiDecodeDate
      public :: cdiEncodeDate
      public :: cdiDecodeTime
      public :: cdiEncodeTime
      public :: streamOpenRead
      public :: streamOpenWrite
      public :: streamOpenAppend
      public :: streamClose
      public :: streamSync
      public :: streamDefVlist
      public :: streamInqVlist
      public :: streamInqFiletype
      public :: streamDefByteorder
      public :: streamInqByteorder
      public :: streamDefZtype
      public :: streamDefZlevel
      public :: streamInqZtype
      public :: streamInqZlevel
      public :: streamDefTimestep
      public :: streamInqTimestep
      public :: streamFilename
      public :: streamFilesuffix
      public :: streamNtsteps
      public :: streamReadVar
      public :: streamWriteVar
      public :: streamReadVarSlice
      public :: streamWriteVarSlice
      public :: streamInqRecord
      public :: streamDefRecord
      public :: streamReadRecord
      public :: streamWriteRecord
      public :: streamCopyRecord
      public :: streamInqGinfo
      public :: vlistCreate
      public :: vlistDestroy
      public :: vlistDuplicate
      public :: vlistCopy
      public :: vlistCopyFlag
      public :: vlistClearFlag
      public :: vlistCat
      public :: vlistMerge
      public :: vlistPrint
      public :: vlistNumber
      public :: vlistNvars
      public :: vlistNgrids
      public :: vlistNzaxis
      public :: vlistDefNtsteps
      public :: vlistNtsteps
      public :: vlistGridsizeMax
      public :: vlistGrid
      public :: vlistGridIndex
      public :: vlistChangeGridIndex
      public :: vlistChangeGrid
      public :: vlistZaxis
      public :: vlistZaxisIndex
      public :: vlistChangeZaxisIndex
      public :: vlistChangeZaxis
      public :: vlistNrecs
      public :: vlistDefTaxis
      public :: vlistInqTaxis
      public :: vlistDefTable
      public :: vlistInqTable
      public :: vlistDefInstitut
      public :: vlistInqInstitut
      public :: vlistDefModel
      public :: vlistInqModel
      public :: vlistDefVar
      public :: vlistChangeVarGrid
      public :: vlistChangeVarZaxis
      public :: vlistInqVar
      public :: vlistInqVarGrid
      public :: vlistInqVarZaxis
      public :: vlistInqVarTime
      public :: vlistDefVarZtype
      public :: vlistInqVarZtype
      public :: vlistDefVarZlevel
      public :: vlistInqVarZlevel
      public :: vlistDefVarParam
      public :: vlistInqVarParam
      public :: vlistDefVarCode
      public :: vlistInqVarCode
      public :: vlistDefVarDatatype
      public :: vlistInqVarDatatype
      public :: vlistInqVarNumber
      public :: vlistDefVarInstitut
      public :: vlistInqVarInstitut
      public :: vlistDefVarModel
      public :: vlistInqVarModel
      public :: vlistDefVarTable
      public :: vlistInqVarTable
      public :: vlistDefVarName
      public :: vlistInqVarName
      public :: vlistDefVarStdname
      public :: vlistInqVarStdname
      public :: vlistDefVarLongname
      public :: vlistInqVarLongname
      public :: vlistDefVarUnits
      public :: vlistInqVarUnits
      public :: vlistDefVarMissval
      public :: vlistInqVarMissval
      public :: vlistDefVarScalefactor
      public :: vlistInqVarScalefactor
      public :: vlistDefVarAddoffset
      public :: vlistInqVarAddoffset
      public :: vlistDefVarTsteptype
      public :: vlistInqVarTsteptype
      public :: vlistDefVarTimave
      public :: vlistInqVarTimave
      public :: vlistDefVarTimaccu
      public :: vlistInqVarTimaccu
      public :: vlistInqVarSize
      public :: vlistInqVarID
      public :: vlistDefIndex
      public :: vlistInqIndex
      public :: vlistDefFlag
      public :: vlistInqFlag
      public :: vlistFindVar
      public :: vlistFindLevel
      public :: vlistMergedVar
      public :: vlistMergedLevel
      public :: vlistInqNatts
      public :: vlistInqAtt
      public :: vlistDelAtt
      public :: vlistDefAttInt
      public :: vlistDefAttFlt
      public :: vlistDefAttTxt
      public :: vlistInqAttInt
      public :: vlistInqAttFlt
      public :: vlistInqAttTxt
      public :: gridName
      public :: gridNamePtr
      public :: gridCompress
      public :: gridDefMask
      public :: gridInqMask
      public :: gridPrint
      public :: gridSize
      public :: gridCreate
      public :: gridDestroy
      public :: gridDuplicate
      public :: gridInqType
      public :: gridInqSize
      public :: gridDefXsize
      public :: gridInqXsize
      public :: gridDefYsize
      public :: gridInqYsize
      public :: gridDefXvals
      public :: gridInqXvals
      public :: gridDefYvals
      public :: gridInqYvals
      public :: gridDefXname
      public :: gridDefXlongname
      public :: gridDefXunits
      public :: gridDefYname
      public :: gridDefYlongname
      public :: gridDefYunits
      public :: gridInqXname
      public :: gridInqXlongname
      public :: gridInqXstdname
      public :: gridInqXunits
      public :: gridInqYname
      public :: gridInqYlongname
      public :: gridInqYstdname
      public :: gridInqYunits
      public :: gridDefPrec
      public :: gridInqPrec
      public :: gridInqXval
      public :: gridInqYval
      public :: gridInqXinc
      public :: gridInqYinc
      public :: gridIsCircular
      public :: gridIsRotated
      public :: gridInqXpole
      public :: gridDefXpole
      public :: gridInqYpole
      public :: gridDefYpole
      public :: gridInqAngle
      public :: gridDefAngle
      public :: gridDefTrunc
      public :: gridInqTrunc
      public :: gridInqGMEnd
      public :: gridDefGMEnd
      public :: gridInqGMEni
      public :: gridDefGMEni
      public :: gridInqGMEni2
      public :: gridDefGMEni2
      public :: gridInqGMEni3
      public :: gridDefGMEni3
      public :: gridDefLCC
      public :: gridInqLCC
      public :: gridDefLcc2
      public :: gridInqLcc2
      public :: gridDefLaea
      public :: gridInqLaea
      public :: gridDefArea
      public :: gridInqArea
      public :: gridHasArea
      public :: gridDefNvertex
      public :: gridInqNvertex
      public :: gridDefXbounds
      public :: gridInqXbounds
      public :: gridDefYbounds
      public :: gridInqYbounds
      public :: gridDefRowlon
      public :: gridInqRowlon
      public :: gridChangeType
      public :: gridDefComplexPacking
      public :: gridInqComplexPacking
      public :: zaxisName
      public :: zaxisCreate
      public :: zaxisDestroy
      public :: zaxisInqType
      public :: zaxisInqSize
      public :: zaxisDuplicate
      public :: zaxisResize
      public :: zaxisPrint
      public :: zaxisSize
      public :: zaxisDefLevels
      public :: zaxisInqLevels
      public :: zaxisDefLevel
      public :: zaxisInqLevel
      public :: zaxisDefName
      public :: zaxisDefLongname
      public :: zaxisDefUnits
      public :: zaxisInqName
      public :: zaxisInqLongname
      public :: zaxisInqUnits
      public :: zaxisDefPrec
      public :: zaxisInqPrec
      public :: zaxisDefLtype
      public :: zaxisInqLtype
      public :: zaxisInqLevelsPtr
      public :: zaxisDefVct
      public :: zaxisInqVctSize
      public :: zaxisInqVctPtr
      public :: zaxisInqLbounds
      public :: zaxisInqUbounds
      public :: zaxisInqWeights
      public :: zaxisInqLbound
      public :: zaxisInqUbound
      public :: zaxisDefLbounds
      public :: zaxisDefUbounds
      public :: zaxisDefWeights
      public :: zaxisChangeType
      public :: taxisCreate
      public :: taxisDestroy
      public :: taxisDuplicate
      public :: taxisCopyTimestep
      public :: taxisDefType
      public :: taxisDefVdate
      public :: taxisDefVtime
      public :: taxisDefRdate
      public :: taxisDefRtime
      public :: taxisHasBounds
      public :: taxisDeleteBounds
      public :: taxisDefVdateBounds
      public :: taxisDefVtimeBounds
      public :: taxisInqVdateBounds
      public :: taxisInqVtimeBounds
      public :: taxisDefCalendar
      public :: taxisDefTunit
      public :: taxisDefNumavg
      public :: taxisInqType
      public :: taxisInqVdate
      public :: taxisInqVtime
      public :: taxisInqRdate
      public :: taxisInqRtime
      public :: taxisInqCalendar
      public :: taxisInqTunit
      public :: taxisInqNumavg
      public :: tunitNamePtr
      public :: institutDef
      public :: institutInq
      public :: institutInqNumber
      public :: institutInqCenter
      public :: institutInqSubcenter
      public :: institutInqNamePtr
      public :: institutInqLongnamePtr
      public :: modelDef
      public :: modelInq
      public :: modelInqInstitut
      public :: modelInqGribID
      public :: modelInqNamePtr
      public :: tableWriteC
      public :: tableWrite
      public :: tableRead
      public :: tableDef
      public :: tableInqNamePtr
      public :: tableDefEntry
      public :: tableInq
      public :: tableInqNumber
      public :: tableInqNum
      public :: tableInqModel
      public :: tableInqPar
      public :: tableInqParCode
      public :: tableInqParName
      public :: tableInqParLongname
      public :: tableInqParUnits
      public :: tableInqParNamePtr
      public :: tableInqParLongnamePtr
      public :: tableInqParUnitsPtr
      public :: streamDefHistory
      public :: streamInqHistorySize
      public :: streamInqHistoryString
      public :: ctrim

      public :: CDI_UNDEFID
      public :: CDI_GLOBAL
      public :: CDI_BIGENDIAN
      public :: CDI_LITTLEENDIAN
      public :: CDI_REAL
      public :: CDI_COMP
      public :: CDI_BOTH
      public :: CDI_ESYSTEM
      public :: CDI_EINVAL
      public :: CDI_EUFTYPE
      public :: CDI_ELIBNAVAIL
      public :: CDI_EUFSTRUCT
      public :: CDI_EUNC4
      public :: CDI_ELIMIT
      public :: FILETYPE_GRB
      public :: FILETYPE_GRB2
      public :: FILETYPE_NC
      public :: FILETYPE_NC2
      public :: FILETYPE_NC4
      public :: FILETYPE_SRV
      public :: FILETYPE_EXT
      public :: FILETYPE_IEG
      public :: COMPRESS_NONE
      public :: COMPRESS_SZIP
      public :: COMPRESS_GZIP
      public :: COMPRESS_BZIP2
      public :: COMPRESS_ZIP
      public :: COMPRESS_JPEG
      public :: DATATYPE_PACK
      public :: DATATYPE_PACK1
      public :: DATATYPE_PACK2
      public :: DATATYPE_PACK3
      public :: DATATYPE_PACK4
      public :: DATATYPE_PACK5
      public :: DATATYPE_PACK6
      public :: DATATYPE_PACK7
      public :: DATATYPE_PACK8
      public :: DATATYPE_PACK9
      public :: DATATYPE_PACK10
      public :: DATATYPE_PACK11
      public :: DATATYPE_PACK12
      public :: DATATYPE_PACK13
      public :: DATATYPE_PACK14
      public :: DATATYPE_PACK15
      public :: DATATYPE_PACK16
      public :: DATATYPE_PACK17
      public :: DATATYPE_PACK18
      public :: DATATYPE_PACK19
      public :: DATATYPE_PACK20
      public :: DATATYPE_PACK21
      public :: DATATYPE_PACK22
      public :: DATATYPE_PACK23
      public :: DATATYPE_PACK24
      public :: DATATYPE_PACK25
      public :: DATATYPE_PACK26
      public :: DATATYPE_PACK27
      public :: DATATYPE_PACK28
      public :: DATATYPE_PACK29
      public :: DATATYPE_PACK30
      public :: DATATYPE_PACK31
      public :: DATATYPE_PACK32
      public :: DATATYPE_CPX32
      public :: DATATYPE_CPX64
      public :: DATATYPE_FLT32
      public :: DATATYPE_FLT64
      public :: DATATYPE_INT8
      public :: DATATYPE_INT16
      public :: DATATYPE_INT32
      public :: DATATYPE_UINT8
      public :: DATATYPE_UINT16
      public :: DATATYPE_UINT32
      public :: DATATYPE_INT
      public :: DATATYPE_FLT
      public :: DATATYPE_TXT
      public :: DATATYPE_CPX
      public :: GRID_GENERIC
      public :: GRID_GAUSSIAN
      public :: GRID_GAUSSIAN_REDUCED
      public :: GRID_LONLAT
      public :: GRID_SPECTRAL
      public :: GRID_FOURIER
      public :: GRID_GME
      public :: GRID_TRAJECTORY
      public :: GRID_UNSTRUCTURED
      public :: GRID_CURVILINEAR
      public :: GRID_LCC
      public :: GRID_LCC2
      public :: GRID_LAEA
      public :: GRID_SINUSOIDAL
      public :: ZAXIS_SURFACE
      public :: ZAXIS_GENERIC
      public :: ZAXIS_HYBRID
      public :: ZAXIS_HYBRID_HALF
      public :: ZAXIS_PRESSURE
      public :: ZAXIS_HEIGHT
      public :: ZAXIS_DEPTH_BELOW_SEA
      public :: ZAXIS_DEPTH_BELOW_LAND
      public :: ZAXIS_ISENTROPIC
      public :: ZAXIS_TRAJECTORY
      public :: ZAXIS_ALTITUDE
      public :: ZAXIS_SIGMA
      public :: ZAXIS_MEANSEA
      public :: TAXIS_ABSOLUTE
      public :: TAXIS_RELATIVE
      public :: TIME_CONSTANT
      public :: TIME_VARIABLE
      public :: TUNIT_SECOND
      public :: TUNIT_MINUTE
      public :: TUNIT_HOUR
      public :: TUNIT_DAY
      public :: TUNIT_MONTH
      public :: TUNIT_YEAR
      public :: TUNIT_QUARTER
      public :: TUNIT_3HOURS
      public :: TUNIT_6HOURS
      public :: TUNIT_12HOURS
      public :: TSTEP_INSTANT
      public :: TSTEP_AVG
      public :: TSTEP_ACCUM
      public :: TSTEP_MAX
      public :: TSTEP_MIN
      public :: TSTEP_DIFF
      public :: TSTEP_RANGE
      public :: TSTEP_INSTANT2
      public :: TSTEP_INSTANT3
      public :: CALENDAR_STANDARD
      public :: CALENDAR_PROLEPTIC
      public :: CALENDAR_360DAYS
      public :: CALENDAR_365DAYS
      public :: CALENDAR_366DAYS
      public :: CALENDAR_NONE

contains

    subroutine ctrim(str)
    use iso_c_binding
    character(kind=c_char, len=*)  :: str
    character                      :: c
    integer                        :: i

    do i=1,len(str)
      c = str(i:i)
      if (c == c_null_char) then
        str(i:len(str)) = ' '
      end if
    end do

    end subroutine ctrim
  

end module mo_cdi
  

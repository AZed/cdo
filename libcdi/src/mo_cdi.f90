
module mo_cdi
      use, intrinsic :: iso_c_binding

      implicit none

      private
  
      integer, parameter :: CDI_MAX_NAME = 256
      integer, parameter :: CDI_UNDEFID = -1
      integer, parameter :: CDI_GLOBAL = -1
      integer, parameter :: CDI_BIGENDIAN = 0
      integer, parameter :: CDI_LITTLEENDIAN = 1
      integer, parameter :: CDI_REAL = 1
      integer, parameter :: CDI_COMP = 2
      integer, parameter :: CDI_BOTH = 3
      integer, parameter :: CDI_ESYSTEM = -10
      integer, parameter :: CDI_EINVAL = -20
      integer, parameter :: CDI_EUFTYPE = -21
      integer, parameter :: CDI_ELIBNAVAIL = -22
      integer, parameter :: CDI_EUFSTRUCT = -23
      integer, parameter :: CDI_EUNC4 = -24
      integer, parameter :: CDI_ELIMIT = -99
      integer, parameter :: FILETYPE_UNDEF = -1
      integer, parameter :: FILETYPE_GRB = 1
      integer, parameter :: FILETYPE_GRB2 = 2
      integer, parameter :: FILETYPE_NC = 3
      integer, parameter :: FILETYPE_NC2 = 4
      integer, parameter :: FILETYPE_NC4 = 5
      integer, parameter :: FILETYPE_NC4C = 6
      integer, parameter :: FILETYPE_SRV = 7
      integer, parameter :: FILETYPE_EXT = 8
      integer, parameter :: FILETYPE_IEG = 9
      integer, parameter :: COMPRESS_NONE = 0
      integer, parameter :: COMPRESS_SZIP = 1
      integer, parameter :: COMPRESS_GZIP = 2
      integer, parameter :: COMPRESS_BZIP2 = 3
      integer, parameter :: COMPRESS_ZIP = 4
      integer, parameter :: COMPRESS_JPEG = 5
      integer, parameter :: DATATYPE_PACK = 0
      integer, parameter :: DATATYPE_PACK1 = 1
      integer, parameter :: DATATYPE_PACK2 = 2
      integer, parameter :: DATATYPE_PACK3 = 3
      integer, parameter :: DATATYPE_PACK4 = 4
      integer, parameter :: DATATYPE_PACK5 = 5
      integer, parameter :: DATATYPE_PACK6 = 6
      integer, parameter :: DATATYPE_PACK7 = 7
      integer, parameter :: DATATYPE_PACK8 = 8
      integer, parameter :: DATATYPE_PACK9 = 9
      integer, parameter :: DATATYPE_PACK10 = 10
      integer, parameter :: DATATYPE_PACK11 = 11
      integer, parameter :: DATATYPE_PACK12 = 12
      integer, parameter :: DATATYPE_PACK13 = 13
      integer, parameter :: DATATYPE_PACK14 = 14
      integer, parameter :: DATATYPE_PACK15 = 15
      integer, parameter :: DATATYPE_PACK16 = 16
      integer, parameter :: DATATYPE_PACK17 = 17
      integer, parameter :: DATATYPE_PACK18 = 18
      integer, parameter :: DATATYPE_PACK19 = 19
      integer, parameter :: DATATYPE_PACK20 = 20
      integer, parameter :: DATATYPE_PACK21 = 21
      integer, parameter :: DATATYPE_PACK22 = 22
      integer, parameter :: DATATYPE_PACK23 = 23
      integer, parameter :: DATATYPE_PACK24 = 24
      integer, parameter :: DATATYPE_PACK25 = 25
      integer, parameter :: DATATYPE_PACK26 = 26
      integer, parameter :: DATATYPE_PACK27 = 27
      integer, parameter :: DATATYPE_PACK28 = 28
      integer, parameter :: DATATYPE_PACK29 = 29
      integer, parameter :: DATATYPE_PACK30 = 30
      integer, parameter :: DATATYPE_PACK31 = 31
      integer, parameter :: DATATYPE_PACK32 = 32
      integer, parameter :: DATATYPE_CPX32 = 64
      integer, parameter :: DATATYPE_CPX64 = 128
      integer, parameter :: DATATYPE_FLT32 = 132
      integer, parameter :: DATATYPE_FLT64 = 164
      integer, parameter :: DATATYPE_INT8 = 208
      integer, parameter :: DATATYPE_INT16 = 216
      integer, parameter :: DATATYPE_INT32 = 232
      integer, parameter :: DATATYPE_UINT8 = 308
      integer, parameter :: DATATYPE_UINT16 = 316
      integer, parameter :: DATATYPE_UINT32 = 332
      integer, parameter :: DATATYPE_INT = 251
      integer, parameter :: DATATYPE_FLT = 252
      integer, parameter :: DATATYPE_TXT = 253
      integer, parameter :: DATATYPE_CPX = 254
      integer, parameter :: DATATYPE_UCHAR = 255
      integer, parameter :: CHUNK_AUTO = 1
      integer, parameter :: CHUNK_GRID = 2
      integer, parameter :: CHUNK_LINES = 3
      integer, parameter :: GRID_GENERIC = 1
      integer, parameter :: GRID_GAUSSIAN = 2
      integer, parameter :: GRID_GAUSSIAN_REDUCED = 3
      integer, parameter :: GRID_LONLAT = 4
      integer, parameter :: GRID_SPECTRAL = 5
      integer, parameter :: GRID_FOURIER = 6
      integer, parameter :: GRID_GME = 7
      integer, parameter :: GRID_TRAJECTORY = 8
      integer, parameter :: GRID_UNSTRUCTURED = 9
      integer, parameter :: GRID_CURVILINEAR = 10
      integer, parameter :: GRID_LCC = 11
      integer, parameter :: GRID_LCC2 = 12
      integer, parameter :: GRID_LAEA = 13
      integer, parameter :: GRID_SINUSOIDAL = 14
      integer, parameter :: GRID_PROJECTION = 15
      integer, parameter :: ZAXIS_SURFACE = 0
      integer, parameter :: ZAXIS_GENERIC = 1
      integer, parameter :: ZAXIS_HYBRID = 2
      integer, parameter :: ZAXIS_HYBRID_HALF = 3
      integer, parameter :: ZAXIS_PRESSURE = 4
      integer, parameter :: ZAXIS_HEIGHT = 5
      integer, parameter :: ZAXIS_DEPTH_BELOW_SEA = 6
      integer, parameter :: ZAXIS_DEPTH_BELOW_LAND = 7
      integer, parameter :: ZAXIS_ISENTROPIC = 8
      integer, parameter :: ZAXIS_TRAJECTORY = 9
      integer, parameter :: ZAXIS_ALTITUDE = 10
      integer, parameter :: ZAXIS_SIGMA = 11
      integer, parameter :: ZAXIS_MEANSEA = 12
      integer, parameter :: ZAXIS_TOA = 13
      integer, parameter :: ZAXIS_SEA_BOTTOM = 14
      integer, parameter :: ZAXIS_ATMOSPHERE = 15
      integer, parameter :: ZAXIS_CLOUD_BASE = 16
      integer, parameter :: ZAXIS_CLOUD_TOP = 17
      integer, parameter :: ZAXIS_ISOTHERM_ZERO = 18
      integer, parameter :: ZAXIS_SNOW = 19
      integer, parameter :: ZAXIS_LAKE_BOTTOM = 20
      integer, parameter :: ZAXIS_SEDIMENT_BOTTOM = 21
      integer, parameter :: ZAXIS_SEDIMENT_BOTTOM_TA = 22
      integer, parameter :: ZAXIS_SEDIMENT_BOTTOM_TW = 23
      integer, parameter :: ZAXIS_MIX_LAYER = 24
      integer, parameter :: ZAXIS_REFERENCE = 25
      integer, parameter :: TIME_CONSTANT = 0
      integer, parameter :: TIME_VARIABLE = 1
      integer, parameter :: TSTEP_CONSTANT = 0
      integer, parameter :: TSTEP_INSTANT = 1
      integer, parameter :: TSTEP_AVG = 2
      integer, parameter :: TSTEP_ACCUM = 3
      integer, parameter :: TSTEP_MAX = 4
      integer, parameter :: TSTEP_MIN = 5
      integer, parameter :: TSTEP_DIFF = 6
      integer, parameter :: TSTEP_RMS = 7
      integer, parameter :: TSTEP_SD = 8
      integer, parameter :: TSTEP_COV = 9
      integer, parameter :: TSTEP_RATIO = 10
      integer, parameter :: TSTEP_RANGE = 11
      integer, parameter :: TSTEP_INSTANT2 = 12
      integer, parameter :: TSTEP_INSTANT3 = 13
      integer, parameter :: TAXIS_ABSOLUTE = 1
      integer, parameter :: TAXIS_RELATIVE = 2
      integer, parameter :: TAXIS_FORECAST = 3
      integer, parameter :: TUNIT_SECOND = 1
      integer, parameter :: TUNIT_MINUTE = 2
      integer, parameter :: TUNIT_QUARTER = 3
      integer, parameter :: TUNIT_30MINUTES = 4
      integer, parameter :: TUNIT_HOUR = 5
      integer, parameter :: TUNIT_3HOURS = 6
      integer, parameter :: TUNIT_6HOURS = 7
      integer, parameter :: TUNIT_12HOURS = 8
      integer, parameter :: TUNIT_DAY = 9
      integer, parameter :: TUNIT_MONTH = 10
      integer, parameter :: TUNIT_YEAR = 11
      integer, parameter :: CALENDAR_STANDARD = 0
      integer, parameter :: CALENDAR_PROLEPTIC = 1
      integer, parameter :: CALENDAR_360DAYS = 2
      integer, parameter :: CALENDAR_365DAYS = 3
      integer, parameter :: CALENDAR_366DAYS = 4
      integer, parameter :: CALENDAR_NONE = 5
      integer, parameter :: CDI_UUID_SIZE = 16
      interface
        function strlen(s) bind(c,name='strlen')
          import :: c_ptr,c_size_t
          type(c_ptr), value :: s
          integer(kind=c_size_t) :: strlen
        end function strlen
      end interface
      interface
        function getchar() bind(c,name='getchar')
          import :: c_int
          integer(kind=c_int) :: getchar
        end function getchar
      end interface
      interface
        function getchar_unlocked() bind(c,name='getchar_unlocked')
          import :: c_int
          integer(kind=c_int) :: getchar_unlocked
        end function getchar_unlocked
      end interface
      interface
        subroutine cdiReset() bind(c,name='cdiReset')
        end subroutine cdiReset
      end interface
      interface
        subroutine cdiDebug(debug) bind(c,name='cdiDebug')
          import :: c_int
          integer(kind=c_int), value :: debug
        end subroutine cdiDebug
      end interface
      interface
        subroutine cdiPrintVersion() bind(c,name='cdiPrintVersion')
        end subroutine cdiPrintVersion
      end interface
      interface
        function cdiHaveFiletype(filetype) bind(c,name='cdiHaveFiletype')
          import :: c_int
          integer(kind=c_int), value :: filetype
          integer(kind=c_int) :: cdiHaveFiletype
        end function cdiHaveFiletype
      end interface
      interface
        subroutine cdiDefMissval(missval) bind(c,name='cdiDefMissval')
          import :: c_double
          real(kind=c_double), value :: missval
        end subroutine cdiDefMissval
      end interface
      interface
        function cdiInqMissval() bind(c,name='cdiInqMissval')
          import :: c_double
          real(kind=c_double) :: cdiInqMissval
        end function cdiInqMissval
      end interface
      interface
        subroutine cdiDefGlobal(string,val) bind(c,name='cdiDefGlobal')
          import :: c_char,c_int
          character(kind=c_char), dimension(*) :: string
          integer(kind=c_int), value :: val
        end subroutine cdiDefGlobal
      end interface
      interface
        function namespaceNew() bind(c,name='namespaceNew')
          import :: c_int
          integer(kind=c_int) :: namespaceNew
        end function namespaceNew
      end interface
      interface
        subroutine namespaceSetActive(namespaceID) bind(c,name='namespaceSetActive')
          import :: c_int
          integer(kind=c_int), value :: namespaceID
        end subroutine namespaceSetActive
      end interface
      interface
        subroutine namespaceDelete(namespaceID) bind(c,name='namespaceDelete')
          import :: c_int
          integer(kind=c_int), value :: namespaceID
        end subroutine namespaceDelete
      end interface
      interface
        subroutine cdiParamToString(param,paramstr,maxlen) bind(c,name='cdiParamToString')
          import :: c_int,c_char
          integer(kind=c_int), value :: param
          character(kind=c_char), dimension(*) :: paramstr
          integer(kind=c_int), value :: maxlen
        end subroutine cdiParamToString
      end interface
      interface
        subroutine cdiDecodeParam(param,pnum,pcat,pdis) bind(c,name='cdiDecodeParam')
          import :: c_int
          integer(kind=c_int), value :: param
          integer(kind=c_int), intent(out) :: pnum
          integer(kind=c_int), intent(out) :: pcat
          integer(kind=c_int), intent(out) :: pdis
        end subroutine cdiDecodeParam
      end interface
      interface
        function cdiEncodeParam(pnum,pcat,pdis) bind(c,name='cdiEncodeParam')
          import :: c_int
          integer(kind=c_int), value :: pnum
          integer(kind=c_int), value :: pcat
          integer(kind=c_int), value :: pdis
          integer(kind=c_int) :: cdiEncodeParam
        end function cdiEncodeParam
      end interface
      interface
        subroutine cdiDecodeDate(date,year,month,day) bind(c,name='cdiDecodeDate')
          import :: c_int
          integer(kind=c_int), value :: date
          integer(kind=c_int), intent(out) :: year
          integer(kind=c_int), intent(out) :: month
          integer(kind=c_int), intent(out) :: day
        end subroutine cdiDecodeDate
      end interface
      interface
        function cdiEncodeDate(year,month,day) bind(c,name='cdiEncodeDate')
          import :: c_int
          integer(kind=c_int), value :: year
          integer(kind=c_int), value :: month
          integer(kind=c_int), value :: day
          integer(kind=c_int) :: cdiEncodeDate
        end function cdiEncodeDate
      end interface
      interface
        subroutine cdiDecodeTime(time,hour,minute,second) bind(c,name='cdiDecodeTime')
          import :: c_int
          integer(kind=c_int), value :: time
          integer(kind=c_int), intent(out) :: hour
          integer(kind=c_int), intent(out) :: minute
          integer(kind=c_int), intent(out) :: second
        end subroutine cdiDecodeTime
      end interface
      interface
        function cdiEncodeTime(hour,minute,second) bind(c,name='cdiEncodeTime')
          import :: c_int
          integer(kind=c_int), value :: hour
          integer(kind=c_int), value :: minute
          integer(kind=c_int), value :: second
          integer(kind=c_int) :: cdiEncodeTime
        end function cdiEncodeTime
      end interface
      interface
        function cdiGetFiletype(path,byteorder) bind(c,name='cdiGetFiletype')
          import :: c_char,c_int
          character(kind=c_char), dimension(*) :: path
          integer(kind=c_int), intent(out) :: byteorder
          integer(kind=c_int) :: cdiGetFiletype
        end function cdiGetFiletype
      end interface
      interface
        function streamOpenRead(path) bind(c,name='streamOpenRead')
          import :: c_char,c_int
          character(kind=c_char), dimension(*) :: path
          integer(kind=c_int) :: streamOpenRead
        end function streamOpenRead
      end interface
      interface
        function streamOpenWrite(path,filetype) bind(c,name='streamOpenWrite')
          import :: c_char,c_int
          character(kind=c_char), dimension(*) :: path
          integer(kind=c_int), value :: filetype
          integer(kind=c_int) :: streamOpenWrite
        end function streamOpenWrite
      end interface
      interface
        function streamOpenAppend(path) bind(c,name='streamOpenAppend')
          import :: c_char,c_int
          character(kind=c_char), dimension(*) :: path
          integer(kind=c_int) :: streamOpenAppend
        end function streamOpenAppend
      end interface
      interface
        subroutine streamClose(streamID) bind(c,name='streamClose')
          import :: c_int
          integer(kind=c_int), value :: streamID
        end subroutine streamClose
      end interface
      interface
        subroutine streamSync(streamID) bind(c,name='streamSync')
          import :: c_int
          integer(kind=c_int), value :: streamID
        end subroutine streamSync
      end interface
      interface
        subroutine streamDefVlist(streamID,vlistID) bind(c,name='streamDefVlist')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: vlistID
        end subroutine streamDefVlist
      end interface
      interface
        function streamInqVlist(streamID) bind(c,name='streamInqVlist')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int) :: streamInqVlist
        end function streamInqVlist
      end interface
      interface
        function streamInqVlistIDorig(streamID) bind(c,name='streamInqVlistIDorig')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int) :: streamInqVlistIDorig
        end function streamInqVlistIDorig
      end interface
      interface
        function streamInqFiletype(streamID) bind(c,name='streamInqFiletype')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int) :: streamInqFiletype
        end function streamInqFiletype
      end interface
      interface
        subroutine streamDefByteorder(streamID,byteorder) bind(c,name='streamDefByteorder')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: byteorder
        end subroutine streamDefByteorder
      end interface
      interface
        function streamInqByteorder(streamID) bind(c,name='streamInqByteorder')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int) :: streamInqByteorder
        end function streamInqByteorder
      end interface
      interface
        subroutine streamDefCompType(streamID,comptype) bind(c,name='streamDefCompType')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: comptype
        end subroutine streamDefCompType
      end interface
      interface
        function streamInqCompType(streamID) bind(c,name='streamInqCompType')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int) :: streamInqCompType
        end function streamInqCompType
      end interface
      interface
        subroutine streamDefCompLevel(streamID,complevel) bind(c,name='streamDefCompLevel')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: complevel
        end subroutine streamDefCompLevel
      end interface
      interface
        function streamInqCompLevel(streamID) bind(c,name='streamInqCompLevel')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int) :: streamInqCompLevel
        end function streamInqCompLevel
      end interface
      interface
        function streamDefTimestep(streamID,tsID) bind(c,name='streamDefTimestep')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: tsID
          integer(kind=c_int) :: streamDefTimestep
        end function streamDefTimestep
      end interface
      interface
        function streamInqTimestep(streamID,tsID) bind(c,name='streamInqTimestep')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: tsID
          integer(kind=c_int) :: streamInqTimestep
        end function streamInqTimestep
      end interface
      interface
        function streamInqCurTimestepID(streamID) bind(c,name='streamInqCurTimestepID')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int) :: streamInqCurTimestepID
        end function streamInqCurTimestepID
      end interface
      interface
        function streamInqNvars(streamID) bind(c,name='streamInqNvars')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int) :: streamInqNvars
        end function streamInqNvars
      end interface
      interface
        subroutine streamWriteVar(streamID,varID,data_vec,nmiss) bind(c,name='streamWriteVar')
          import :: c_int,c_double
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: varID
          real(kind=c_double), intent(in),dimension(*) :: data_vec
          integer(kind=c_int), value :: nmiss
        end subroutine streamWriteVar
      end interface
      interface
        subroutine streamWriteVarF(streamID,varID,data_vec,nmiss) bind(c,name='streamWriteVarF')
          import :: c_int,c_float
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: varID
          real(kind=c_float), intent(in),dimension(*) :: data_vec
          integer(kind=c_int), value :: nmiss
        end subroutine streamWriteVarF
      end interface
      interface
        subroutine streamReadVar(streamID,varID,data_vec,nmiss) bind(c,name='streamReadVar')
          import :: c_int,c_double
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: varID
          real(kind=c_double), intent(out),dimension(*) :: data_vec
          integer(kind=c_int), intent(out) :: nmiss
        end subroutine streamReadVar
      end interface
      interface
        subroutine streamReadVarF(streamID,varID,data_vec,nmiss) bind(c,name='streamReadVarF')
          import :: c_int,c_float
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: varID
          real(kind=c_float), intent(out),dimension(*) :: data_vec
          integer(kind=c_int), intent(out) :: nmiss
        end subroutine streamReadVarF
      end interface
      interface
        subroutine streamWriteVarSlice(streamID,varID,levelID,data_vec,nmiss) bind(c,name='streamWriteVarSlice')
          import :: c_int,c_double
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: levelID
          real(kind=c_double), intent(in),dimension(*) :: data_vec
          integer(kind=c_int), value :: nmiss
        end subroutine streamWriteVarSlice
      end interface
      interface
        subroutine streamWriteVarSliceF(streamID,varID,levelID,data_vec,nmiss) bind(c,name='streamWriteVarSliceF')
          import :: c_int,c_float
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: levelID
          real(kind=c_float), intent(in),dimension(*) :: data_vec
          integer(kind=c_int), value :: nmiss
        end subroutine streamWriteVarSliceF
      end interface
      interface
        subroutine streamReadVarSlice(streamID,varID,levelID,data_vec,nmiss) bind(c,name='streamReadVarSlice')
          import :: c_int,c_double
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: levelID
          real(kind=c_double), intent(out),dimension(*) :: data_vec
          integer(kind=c_int), intent(out) :: nmiss
        end subroutine streamReadVarSlice
      end interface
      interface
        subroutine streamReadVarSliceF(streamID,varID,levelID,data_vec,nmiss) bind(c,name='streamReadVarSliceF')
          import :: c_int,c_float
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: levelID
          real(kind=c_float), intent(out),dimension(*) :: data_vec
          integer(kind=c_int), intent(out) :: nmiss
        end subroutine streamReadVarSliceF
      end interface
      interface
        subroutine streamDefRecord(streamID,varID,levelID) bind(c,name='streamDefRecord')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: levelID
        end subroutine streamDefRecord
      end interface
      interface
        subroutine streamInqRecord(streamID,varID,levelID) bind(c,name='streamInqRecord')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), intent(out) :: varID
          integer(kind=c_int), intent(out) :: levelID
        end subroutine streamInqRecord
      end interface
      interface
        subroutine streamWriteRecord(streamID,data_vec,nmiss) bind(c,name='streamWriteRecord')
          import :: c_int,c_double
          integer(kind=c_int), value :: streamID
          real(kind=c_double), intent(in),dimension(*) :: data_vec
          integer(kind=c_int), value :: nmiss
        end subroutine streamWriteRecord
      end interface
      interface
        subroutine streamWriteRecordF(streamID,data_vec,nmiss) bind(c,name='streamWriteRecordF')
          import :: c_int,c_float
          integer(kind=c_int), value :: streamID
          real(kind=c_float), intent(in),dimension(*) :: data_vec
          integer(kind=c_int), value :: nmiss
        end subroutine streamWriteRecordF
      end interface
      interface
        subroutine streamReadRecord(streamID,data_vec,nmiss) bind(c,name='streamReadRecord')
          import :: c_int,c_double
          integer(kind=c_int), value :: streamID
          real(kind=c_double), intent(out),dimension(*) :: data_vec
          integer(kind=c_int), intent(out) :: nmiss
        end subroutine streamReadRecord
      end interface
      interface
        subroutine streamCopyRecord(streamIDdest,streamIDsrc) bind(c,name='streamCopyRecord')
          import :: c_int
          integer(kind=c_int), value :: streamIDdest
          integer(kind=c_int), value :: streamIDsrc
        end subroutine streamCopyRecord
      end interface
      interface
        function vlistCreate() bind(c,name='vlistCreate')
          import :: c_int
          integer(kind=c_int) :: vlistCreate
        end function vlistCreate
      end interface
      interface
        subroutine vlistDestroy(vlistID) bind(c,name='vlistDestroy')
          import :: c_int
          integer(kind=c_int), value :: vlistID
        end subroutine vlistDestroy
      end interface
      interface
        function vlistDuplicate(vlistID) bind(c,name='vlistDuplicate')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int) :: vlistDuplicate
        end function vlistDuplicate
      end interface
      interface
        subroutine vlistCopy(vlistID2,vlistID1) bind(c,name='vlistCopy')
          import :: c_int
          integer(kind=c_int), value :: vlistID2
          integer(kind=c_int), value :: vlistID1
        end subroutine vlistCopy
      end interface
      interface
        subroutine vlistCopyFlag(vlistID2,vlistID1) bind(c,name='vlistCopyFlag')
          import :: c_int
          integer(kind=c_int), value :: vlistID2
          integer(kind=c_int), value :: vlistID1
        end subroutine vlistCopyFlag
      end interface
      interface
        subroutine vlistClearFlag(vlistID) bind(c,name='vlistClearFlag')
          import :: c_int
          integer(kind=c_int), value :: vlistID
        end subroutine vlistClearFlag
      end interface
      interface
        subroutine vlistCat(vlistID2,vlistID1) bind(c,name='vlistCat')
          import :: c_int
          integer(kind=c_int), value :: vlistID2
          integer(kind=c_int), value :: vlistID1
        end subroutine vlistCat
      end interface
      interface
        subroutine vlistMerge(vlistID2,vlistID1) bind(c,name='vlistMerge')
          import :: c_int
          integer(kind=c_int), value :: vlistID2
          integer(kind=c_int), value :: vlistID1
        end subroutine vlistMerge
      end interface
      interface
        subroutine vlistPrint(vlistID) bind(c,name='vlistPrint')
          import :: c_int
          integer(kind=c_int), value :: vlistID
        end subroutine vlistPrint
      end interface
      interface
        function vlistNumber(vlistID) bind(c,name='vlistNumber')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int) :: vlistNumber
        end function vlistNumber
      end interface
      interface
        function vlistNvars(vlistID) bind(c,name='vlistNvars')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int) :: vlistNvars
        end function vlistNvars
      end interface
      interface
        function vlistNgrids(vlistID) bind(c,name='vlistNgrids')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int) :: vlistNgrids
        end function vlistNgrids
      end interface
      interface
        function vlistNzaxis(vlistID) bind(c,name='vlistNzaxis')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int) :: vlistNzaxis
        end function vlistNzaxis
      end interface
      interface
        subroutine vlistDefNtsteps(vlistID,nts) bind(c,name='vlistDefNtsteps')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: nts
        end subroutine vlistDefNtsteps
      end interface
      interface
        function vlistNtsteps(vlistID) bind(c,name='vlistNtsteps')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int) :: vlistNtsteps
        end function vlistNtsteps
      end interface
      interface
        function vlistGridsizeMax(vlistID) bind(c,name='vlistGridsizeMax')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int) :: vlistGridsizeMax
        end function vlistGridsizeMax
      end interface
      interface
        function vlistGrid(vlistID,index) bind(c,name='vlistGrid')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: index
          integer(kind=c_int) :: vlistGrid
        end function vlistGrid
      end interface
      interface
        function vlistGridIndex(vlistID,gridID) bind(c,name='vlistGridIndex')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: vlistGridIndex
        end function vlistGridIndex
      end interface
      interface
        subroutine vlistChangeGridIndex(vlistID,index,gridID) bind(c,name='vlistChangeGridIndex')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: index
          integer(kind=c_int), value :: gridID
        end subroutine vlistChangeGridIndex
      end interface
      interface
        subroutine vlistChangeGrid(vlistID,gridID1,gridID2) bind(c,name='vlistChangeGrid')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: gridID1
          integer(kind=c_int), value :: gridID2
        end subroutine vlistChangeGrid
      end interface
      interface
        function vlistZaxis(vlistID,index) bind(c,name='vlistZaxis')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: index
          integer(kind=c_int) :: vlistZaxis
        end function vlistZaxis
      end interface
      interface
        function vlistZaxisIndex(vlistID,zaxisID) bind(c,name='vlistZaxisIndex')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int) :: vlistZaxisIndex
        end function vlistZaxisIndex
      end interface
      interface
        subroutine vlistChangeZaxisIndex(vlistID,index,zaxisID) bind(c,name='vlistChangeZaxisIndex')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: index
          integer(kind=c_int), value :: zaxisID
        end subroutine vlistChangeZaxisIndex
      end interface
      interface
        subroutine vlistChangeZaxis(vlistID,zaxisID1,zaxisID2) bind(c,name='vlistChangeZaxis')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: zaxisID1
          integer(kind=c_int), value :: zaxisID2
        end subroutine vlistChangeZaxis
      end interface
      interface
        function vlistNrecs(vlistID) bind(c,name='vlistNrecs')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int) :: vlistNrecs
        end function vlistNrecs
      end interface
      interface
        subroutine vlistDefTaxis(vlistID,taxisID) bind(c,name='vlistDefTaxis')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: taxisID
        end subroutine vlistDefTaxis
      end interface
      interface
        function vlistInqTaxis(vlistID) bind(c,name='vlistInqTaxis')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int) :: vlistInqTaxis
        end function vlistInqTaxis
      end interface
      interface
        subroutine vlistDefTable(vlistID,tableID) bind(c,name='vlistDefTable')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: tableID
        end subroutine vlistDefTable
      end interface
      interface
        function vlistInqTable(vlistID) bind(c,name='vlistInqTable')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int) :: vlistInqTable
        end function vlistInqTable
      end interface
      interface
        subroutine vlistDefInstitut(vlistID,instID) bind(c,name='vlistDefInstitut')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: instID
        end subroutine vlistDefInstitut
      end interface
      interface
        function vlistInqInstitut(vlistID) bind(c,name='vlistInqInstitut')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int) :: vlistInqInstitut
        end function vlistInqInstitut
      end interface
      interface
        subroutine vlistDefModel(vlistID,modelID) bind(c,name='vlistDefModel')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: modelID
        end subroutine vlistDefModel
      end interface
      interface
        function vlistInqModel(vlistID) bind(c,name='vlistInqModel')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int) :: vlistInqModel
        end function vlistInqModel
      end interface
      interface
        function vlistDefVar(vlistID,gridID,zaxisID,tsteptype) bind(c,name='vlistDefVar')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int), value :: tsteptype
          integer(kind=c_int) :: vlistDefVar
        end function vlistDefVar
      end interface
      interface
        subroutine vlistChangeVarGrid(vlistID,varID,gridID) bind(c,name='vlistChangeVarGrid')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: gridID
        end subroutine vlistChangeVarGrid
      end interface
      interface
        subroutine vlistChangeVarZaxis(vlistID,varID,zaxisID) bind(c,name='vlistChangeVarZaxis')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: zaxisID
        end subroutine vlistChangeVarZaxis
      end interface
      interface
        subroutine vlistInqVar(vlistID,varID,gridID,zaxisID,tsteptype) bind(c,name='vlistInqVar')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), intent(out) :: gridID
          integer(kind=c_int), intent(out) :: zaxisID
          integer(kind=c_int), intent(out) :: tsteptype
        end subroutine vlistInqVar
      end interface
      interface
        function vlistInqVarGrid(vlistID,varID) bind(c,name='vlistInqVarGrid')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarGrid
        end function vlistInqVarGrid
      end interface
      interface
        function vlistInqVarZaxis(vlistID,varID) bind(c,name='vlistInqVarZaxis')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarZaxis
        end function vlistInqVarZaxis
      end interface
      interface
        function vlistInqVarID(vlistID,code) bind(c,name='vlistInqVarID')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: code
          integer(kind=c_int) :: vlistInqVarID
        end function vlistInqVarID
      end interface
      interface
        subroutine vlistDefVarTsteptype(vlistID,varID,tsteptype) bind(c,name='vlistDefVarTsteptype')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: tsteptype
        end subroutine vlistDefVarTsteptype
      end interface
      interface
        function vlistInqVarTsteptype(vlistID,varID) bind(c,name='vlistInqVarTsteptype')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarTsteptype
        end function vlistInqVarTsteptype
      end interface
      interface
        subroutine vlistDefVarCompType(vlistID,varID,comptype) bind(c,name='vlistDefVarCompType')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: comptype
        end subroutine vlistDefVarCompType
      end interface
      interface
        function vlistInqVarCompType(vlistID,varID) bind(c,name='vlistInqVarCompType')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarCompType
        end function vlistInqVarCompType
      end interface
      interface
        subroutine vlistDefVarCompLevel(vlistID,varID,complevel) bind(c,name='vlistDefVarCompLevel')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: complevel
        end subroutine vlistDefVarCompLevel
      end interface
      interface
        function vlistInqVarCompLevel(vlistID,varID) bind(c,name='vlistInqVarCompLevel')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarCompLevel
        end function vlistInqVarCompLevel
      end interface
      interface
        subroutine vlistDefVarParam(vlistID,varID,param) bind(c,name='vlistDefVarParam')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: param
        end subroutine vlistDefVarParam
      end interface
      interface
        function vlistInqVarParam(vlistID,varID) bind(c,name='vlistInqVarParam')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarParam
        end function vlistInqVarParam
      end interface
      interface
        subroutine vlistDefVarCode(vlistID,varID,code) bind(c,name='vlistDefVarCode')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: code
        end subroutine vlistDefVarCode
      end interface
      interface
        function vlistInqVarCode(vlistID,varID) bind(c,name='vlistInqVarCode')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarCode
        end function vlistInqVarCode
      end interface
      interface
        subroutine vlistDefVarDatatype(vlistID,varID,datatype) bind(c,name='vlistDefVarDatatype')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: datatype
        end subroutine vlistDefVarDatatype
      end interface
      interface
        function vlistInqVarDatatype(vlistID,varID) bind(c,name='vlistInqVarDatatype')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarDatatype
        end function vlistInqVarDatatype
      end interface
      interface
        subroutine vlistDefVarChunkType(vlistID,varID,chunktype) bind(c,name='vlistDefVarChunkType')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: chunktype
        end subroutine vlistDefVarChunkType
      end interface
      interface
        function vlistInqVarChunkType(vlistID,varID) bind(c,name='vlistInqVarChunkType')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarChunkType
        end function vlistInqVarChunkType
      end interface
      interface
        subroutine vlistDefVarXYZ(vlistID,varID,xyz) bind(c,name='vlistDefVarXYZ')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: xyz
        end subroutine vlistDefVarXYZ
      end interface
      interface
        function vlistInqVarXYZ(vlistID,varID) bind(c,name='vlistInqVarXYZ')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarXYZ
        end function vlistInqVarXYZ
      end interface
      interface
        function vlistInqVarNumber(vlistID,varID) bind(c,name='vlistInqVarNumber')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarNumber
        end function vlistInqVarNumber
      end interface
      interface
        subroutine vlistDefVarInstitut(vlistID,varID,instID) bind(c,name='vlistDefVarInstitut')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: instID
        end subroutine vlistDefVarInstitut
      end interface
      interface
        function vlistInqVarInstitut(vlistID,varID) bind(c,name='vlistInqVarInstitut')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarInstitut
        end function vlistInqVarInstitut
      end interface
      interface
        subroutine vlistDefVarModel(vlistID,varID,modelID) bind(c,name='vlistDefVarModel')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: modelID
        end subroutine vlistDefVarModel
      end interface
      interface
        function vlistInqVarModel(vlistID,varID) bind(c,name='vlistInqVarModel')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarModel
        end function vlistInqVarModel
      end interface
      interface
        subroutine vlistDefVarTable(vlistID,varID,tableID) bind(c,name='vlistDefVarTable')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: tableID
        end subroutine vlistDefVarTable
      end interface
      interface
        function vlistInqVarTable(vlistID,varID) bind(c,name='vlistInqVarTable')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarTable
        end function vlistInqVarTable
      end interface
      interface
        subroutine vlistDefVarName(vlistID,varID,name) bind(c,name='vlistDefVarName')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: name
        end subroutine vlistDefVarName
      end interface
      interface
        subroutine vlistInqVarName(vlistID,varID,name) bind(c,name='vlistInqVarName')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: name
        end subroutine vlistInqVarName
      end interface
      interface
        subroutine vlistDefVarStdname(vlistID,varID,stdname) bind(c,name='vlistDefVarStdname')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: stdname
        end subroutine vlistDefVarStdname
      end interface
      interface
        subroutine vlistInqVarStdname(vlistID,varID,stdname) bind(c,name='vlistInqVarStdname')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: stdname
        end subroutine vlistInqVarStdname
      end interface
      interface
        subroutine vlistDefVarLongname(vlistID,varID,longname) bind(c,name='vlistDefVarLongname')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: longname
        end subroutine vlistDefVarLongname
      end interface
      interface
        subroutine vlistInqVarLongname(vlistID,varID,longname) bind(c,name='vlistInqVarLongname')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: longname
        end subroutine vlistInqVarLongname
      end interface
      interface
        subroutine vlistDefVarUnits(vlistID,varID,units) bind(c,name='vlistDefVarUnits')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: units
        end subroutine vlistDefVarUnits
      end interface
      interface
        subroutine vlistInqVarUnits(vlistID,varID,units) bind(c,name='vlistInqVarUnits')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: units
        end subroutine vlistInqVarUnits
      end interface
      interface
        subroutine vlistDefVarMissval(vlistID,varID,missval) bind(c,name='vlistDefVarMissval')
          import :: c_int,c_double
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          real(kind=c_double), value :: missval
        end subroutine vlistDefVarMissval
      end interface
      interface
        function vlistInqVarMissval(vlistID,varID) bind(c,name='vlistInqVarMissval')
          import :: c_int,c_double
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          real(kind=c_double) :: vlistInqVarMissval
        end function vlistInqVarMissval
      end interface
      interface
        subroutine vlistDefVarExtra(vlistID,varID,extra) bind(c,name='vlistDefVarExtra')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: extra
        end subroutine vlistDefVarExtra
      end interface
      interface
        subroutine vlistInqVarExtra(vlistID,varID,extra) bind(c,name='vlistInqVarExtra')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: extra
        end subroutine vlistInqVarExtra
      end interface
      interface
        subroutine vlistDefVarScalefactor(vlistID,varID,scalefactor) bind(c,name='vlistDefVarScalefactor')
          import :: c_int,c_double
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          real(kind=c_double), value :: scalefactor
        end subroutine vlistDefVarScalefactor
      end interface
      interface
        function vlistInqVarScalefactor(vlistID,varID) bind(c,name='vlistInqVarScalefactor')
          import :: c_int,c_double
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          real(kind=c_double) :: vlistInqVarScalefactor
        end function vlistInqVarScalefactor
      end interface
      interface
        subroutine vlistDefVarAddoffset(vlistID,varID,addoffset) bind(c,name='vlistDefVarAddoffset')
          import :: c_int,c_double
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          real(kind=c_double), value :: addoffset
        end subroutine vlistDefVarAddoffset
      end interface
      interface
        function vlistInqVarAddoffset(vlistID,varID) bind(c,name='vlistInqVarAddoffset')
          import :: c_int,c_double
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          real(kind=c_double) :: vlistInqVarAddoffset
        end function vlistInqVarAddoffset
      end interface
      interface
        subroutine vlistDefVarTimave(vlistID,varID,timave) bind(c,name='vlistDefVarTimave')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: timave
        end subroutine vlistDefVarTimave
      end interface
      interface
        function vlistInqVarTimave(vlistID,varID) bind(c,name='vlistInqVarTimave')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarTimave
        end function vlistInqVarTimave
      end interface
      interface
        subroutine vlistDefVarTimaccu(vlistID,varID,timaccu) bind(c,name='vlistDefVarTimaccu')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: timaccu
        end subroutine vlistDefVarTimaccu
      end interface
      interface
        function vlistInqVarTimaccu(vlistID,varID) bind(c,name='vlistInqVarTimaccu')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarTimaccu
        end function vlistInqVarTimaccu
      end interface
      interface
        subroutine vlistDefVarTypeOfGeneratingProcess(vlistID,varID,typeOfGeneratingProcess) bind(c,&
       name='vlistDefVarTypeOfGeneratingProcess')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: typeOfGeneratingProcess
        end subroutine vlistDefVarTypeOfGeneratingProcess
      end interface
      interface
        function vlistInqVarTypeOfGeneratingProcess(vlistID,varID) bind(c,name='vlistInqVarTypeOfGeneratingProcess')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarTypeOfGeneratingProcess
        end function vlistInqVarTypeOfGeneratingProcess
      end interface
      interface
        subroutine vlistDefVarProductDefinitionTemplate(vlistID,varID,productDefinitionTemplate) bind(c,&
       name='vlistDefVarProductDefinitionTemplate')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: productDefinitionTemplate
        end subroutine vlistDefVarProductDefinitionTemplate
      end interface
      interface
        function vlistInqVarProductDefinitionTemplate(vlistID,varID) bind(c,name='vlistInqVarProductDefinitionTemplate')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarProductDefinitionTemplate
        end function vlistInqVarProductDefinitionTemplate
      end interface
      interface
        function vlistInqVarSize(vlistID,varID) bind(c,name='vlistInqVarSize')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistInqVarSize
        end function vlistInqVarSize
      end interface
      interface
        subroutine vlistDefIndex(vlistID,varID,levID,index) bind(c,name='vlistDefIndex')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: levID
          integer(kind=c_int), value :: index
        end subroutine vlistDefIndex
      end interface
      interface
        function vlistInqIndex(vlistID,varID,levID) bind(c,name='vlistInqIndex')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: levID
          integer(kind=c_int) :: vlistInqIndex
        end function vlistInqIndex
      end interface
      interface
        subroutine vlistDefFlag(vlistID,varID,levID,flag) bind(c,name='vlistDefFlag')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: levID
          integer(kind=c_int), value :: flag
        end subroutine vlistDefFlag
      end interface
      interface
        function vlistInqFlag(vlistID,varID,levID) bind(c,name='vlistInqFlag')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: levID
          integer(kind=c_int) :: vlistInqFlag
        end function vlistInqFlag
      end interface
      interface
        function vlistFindVar(vlistID,fvarID) bind(c,name='vlistFindVar')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: fvarID
          integer(kind=c_int) :: vlistFindVar
        end function vlistFindVar
      end interface
      interface
        function vlistFindLevel(vlistID,fvarID,flevelID) bind(c,name='vlistFindLevel')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: fvarID
          integer(kind=c_int), value :: flevelID
          integer(kind=c_int) :: vlistFindLevel
        end function vlistFindLevel
      end interface
      interface
        function vlistMergedVar(vlistID,varID) bind(c,name='vlistMergedVar')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int) :: vlistMergedVar
        end function vlistMergedVar
      end interface
      interface
        function vlistMergedLevel(vlistID,varID,levelID) bind(c,name='vlistMergedLevel')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: levelID
          integer(kind=c_int) :: vlistMergedLevel
        end function vlistMergedLevel
      end interface
      interface
        subroutine vlistDefVarEnsemble(vlistID,varID,ensID,ensCount,forecast_type) bind(c,name='vlistDefVarEnsemble')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: ensID
          integer(kind=c_int), value :: ensCount
          integer(kind=c_int), value :: forecast_type
        end subroutine vlistDefVarEnsemble
      end interface
      interface
        function vlistInqVarEnsemble(vlistID,varID,ensID,ensCount,forecast_type) bind(c,name='vlistInqVarEnsemble')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), intent(out) :: ensID
          integer(kind=c_int), intent(out) :: ensCount
          integer(kind=c_int), intent(out) :: forecast_type
          integer(kind=c_int) :: vlistInqVarEnsemble
        end function vlistInqVarEnsemble
      end interface
      interface
        subroutine cdiClearAdditionalKeys() bind(c,name='cdiClearAdditionalKeys')
        end subroutine cdiClearAdditionalKeys
      end interface
      interface
        subroutine cdiDefAdditionalKey(string) bind(c,name='cdiDefAdditionalKey')
          import :: c_char
          character(kind=c_char), dimension(*) :: string
        end subroutine cdiDefAdditionalKey
      end interface
      interface
        subroutine vlistDefVarIntKey(vlistID,varID,name,value) bind(c,name='vlistDefVarIntKey')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int), value :: value
        end subroutine vlistDefVarIntKey
      end interface
      interface
        subroutine vlistDefVarDblKey(vlistID,varID,name,value) bind(c,name='vlistDefVarDblKey')
          import :: c_int,c_char,c_double
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: name
          real(kind=c_double), value :: value
        end subroutine vlistDefVarDblKey
      end interface
      interface
        function vlistHasVarKey(vlistID,varID,name) bind(c,name='vlistHasVarKey')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int) :: vlistHasVarKey
        end function vlistHasVarKey
      end interface
      interface
        function vlistInqVarDblKey(vlistID,varID,name) bind(c,name='vlistInqVarDblKey')
          import :: c_int,c_char,c_double
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: name
          real(kind=c_double) :: vlistInqVarDblKey
        end function vlistInqVarDblKey
      end interface
      interface
        function vlistInqVarIntKey(vlistID,varID,name) bind(c,name='vlistInqVarIntKey')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int) :: vlistInqVarIntKey
        end function vlistInqVarIntKey
      end interface
      interface
        function vlistInqNatts(vlistID,varID,nattsp) bind(c,name='vlistInqNatts')
          import :: c_int
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), intent(out) :: nattsp
          integer(kind=c_int) :: vlistInqNatts
        end function vlistInqNatts
      end interface
      interface
        function vlistInqAtt(vlistID,varID,attrnum,name,typep,lenp) bind(c,name='vlistInqAtt')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          integer(kind=c_int), value :: attrnum
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int), intent(out) :: typep
          integer(kind=c_int), intent(out) :: lenp
          integer(kind=c_int) :: vlistInqAtt
        end function vlistInqAtt
      end interface
      interface
        function vlistDelAtt(vlistID,varID,name) bind(c,name='vlistDelAtt')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int) :: vlistDelAtt
        end function vlistDelAtt
      end interface
      interface
        function vlistDefAttInt(vlistID,varID,name,type,len,ip_vec) bind(c,name='vlistDefAttInt')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int), value :: type
          integer(kind=c_int), value :: len
          integer(kind=c_int), intent(in),dimension(*) :: ip_vec
          integer(kind=c_int) :: vlistDefAttInt
        end function vlistDefAttInt
      end interface
      interface
        function vlistDefAttFlt(vlistID,varID,name,type,len,dp_vec) bind(c,name='vlistDefAttFlt')
          import :: c_int,c_char,c_double
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int), value :: type
          integer(kind=c_int), value :: len
          real(kind=c_double), intent(in),dimension(*) :: dp_vec
          integer(kind=c_int) :: vlistDefAttFlt
        end function vlistDefAttFlt
      end interface
      interface
        function vlistDefAttTxt(vlistID,varID,name,len,tp_cbuf) bind(c,name='vlistDefAttTxt')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int), value :: len
          character(kind=c_char), dimension(*) :: tp_cbuf
          integer(kind=c_int) :: vlistDefAttTxt
        end function vlistDefAttTxt
      end interface
      interface
        function vlistInqAttInt(vlistID,varID,name,mlen,ip_vec) bind(c,name='vlistInqAttInt')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int), value :: mlen
          integer(kind=c_int), intent(out),dimension(*) :: ip_vec
          integer(kind=c_int) :: vlistInqAttInt
        end function vlistInqAttInt
      end interface
      interface
        function vlistInqAttFlt(vlistID,varID,name,mlen,dp_vec) bind(c,name='vlistInqAttFlt')
          import :: c_int,c_char,c_double
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int), value :: mlen
          real(kind=c_double), intent(out),dimension(*) :: dp_vec
          integer(kind=c_int) :: vlistInqAttFlt
        end function vlistInqAttFlt
      end interface
      interface
        function vlistInqAttTxt(vlistID,varID,name,mlen,tp_cbuf) bind(c,name='vlistInqAttTxt')
          import :: c_int,c_char
          integer(kind=c_int), value :: vlistID
          integer(kind=c_int), value :: varID
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int), value :: mlen
          character(kind=c_char), dimension(*) :: tp_cbuf
          integer(kind=c_int) :: vlistInqAttTxt
        end function vlistInqAttTxt
      end interface
      interface
        subroutine gridName(gridtype,gridnamev) bind(c,name='gridName')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridtype
          character(kind=c_char), dimension(*) :: gridnamev
        end subroutine gridName
      end interface
      interface
        subroutine gridCompress(gridID) bind(c,name='gridCompress')
          import :: c_int
          integer(kind=c_int), value :: gridID
        end subroutine gridCompress
      end interface
      interface
        subroutine gridDefMaskGME(gridID,mask_vec) bind(c,name='gridDefMaskGME')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), intent(in),dimension(*) :: mask_vec
        end subroutine gridDefMaskGME
      end interface
      interface
        function gridInqMaskGME(gridID,mask_vec) bind(c,name='gridInqMaskGME')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), intent(out),dimension(*) :: mask_vec
          integer(kind=c_int) :: gridInqMaskGME
        end function gridInqMaskGME
      end interface
      interface
        subroutine gridDefMask(gridID,mask_vec) bind(c,name='gridDefMask')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), intent(in),dimension(*) :: mask_vec
        end subroutine gridDefMask
      end interface
      interface
        function gridInqMask(gridID,mask_vec) bind(c,name='gridInqMask')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), intent(out),dimension(*) :: mask_vec
          integer(kind=c_int) :: gridInqMask
        end function gridInqMask
      end interface
      interface
        subroutine gridPrint(gridID,opt) bind(c,name='gridPrint')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: opt
        end subroutine gridPrint
      end interface
      interface
        function gridCreate(gridtype,size) bind(c,name='gridCreate')
          import :: c_int
          integer(kind=c_int), value :: gridtype
          integer(kind=c_int), value :: size
          integer(kind=c_int) :: gridCreate
        end function gridCreate
      end interface
      interface
        subroutine gridDestroy(gridID) bind(c,name='gridDestroy')
          import :: c_int
          integer(kind=c_int), value :: gridID
        end subroutine gridDestroy
      end interface
      interface
        function gridDuplicate(gridID) bind(c,name='gridDuplicate')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridDuplicate
        end function gridDuplicate
      end interface
      interface
        function gridInqType(gridID) bind(c,name='gridInqType')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqType
        end function gridInqType
      end interface
      interface
        function gridInqSize(gridID) bind(c,name='gridInqSize')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqSize
        end function gridInqSize
      end interface
      interface
        subroutine gridDefXsize(gridID,xsize) bind(c,name='gridDefXsize')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: xsize
        end subroutine gridDefXsize
      end interface
      interface
        function gridInqXsize(gridID) bind(c,name='gridInqXsize')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqXsize
        end function gridInqXsize
      end interface
      interface
        subroutine gridDefYsize(gridID,ysize) bind(c,name='gridDefYsize')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: ysize
        end subroutine gridDefYsize
      end interface
      interface
        function gridInqYsize(gridID) bind(c,name='gridInqYsize')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqYsize
        end function gridInqYsize
      end interface
      interface
        subroutine gridDefNP(gridID,np) bind(c,name='gridDefNP')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: np
        end subroutine gridDefNP
      end interface
      interface
        function gridInqNP(gridID) bind(c,name='gridInqNP')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqNP
        end function gridInqNP
      end interface
      interface
        subroutine gridDefXvals(gridID,xvals_vec) bind(c,name='gridDefXvals')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), intent(in),dimension(*) :: xvals_vec
        end subroutine gridDefXvals
      end interface
      interface
        function gridInqXvals(gridID,xvals_vec) bind(c,name='gridInqXvals')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), intent(out),dimension(*) :: xvals_vec
          integer(kind=c_int) :: gridInqXvals
        end function gridInqXvals
      end interface
      interface
        subroutine gridDefYvals(gridID,yvals_vec) bind(c,name='gridDefYvals')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), intent(in),dimension(*) :: yvals_vec
        end subroutine gridDefYvals
      end interface
      interface
        function gridInqYvals(gridID,yvals_vec) bind(c,name='gridInqYvals')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), intent(out),dimension(*) :: yvals_vec
          integer(kind=c_int) :: gridInqYvals
        end function gridInqYvals
      end interface
      interface
        subroutine gridDefXname(gridID,xname) bind(c,name='gridDefXname')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: xname
        end subroutine gridDefXname
      end interface
      interface
        subroutine gridInqXname(gridID,xname) bind(c,name='gridInqXname')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: xname
        end subroutine gridInqXname
      end interface
      interface
        subroutine gridDefXlongname(gridID,xlongname) bind(c,name='gridDefXlongname')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: xlongname
        end subroutine gridDefXlongname
      end interface
      interface
        subroutine gridInqXlongname(gridID,xlongname) bind(c,name='gridInqXlongname')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: xlongname
        end subroutine gridInqXlongname
      end interface
      interface
        subroutine gridDefXunits(gridID,xunits) bind(c,name='gridDefXunits')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: xunits
        end subroutine gridDefXunits
      end interface
      interface
        subroutine gridInqXunits(gridID,xunits) bind(c,name='gridInqXunits')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: xunits
        end subroutine gridInqXunits
      end interface
      interface
        subroutine gridDefYname(gridID,yname) bind(c,name='gridDefYname')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: yname
        end subroutine gridDefYname
      end interface
      interface
        subroutine gridInqYname(gridID,yname) bind(c,name='gridInqYname')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: yname
        end subroutine gridInqYname
      end interface
      interface
        subroutine gridDefYlongname(gridID,ylongname) bind(c,name='gridDefYlongname')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: ylongname
        end subroutine gridDefYlongname
      end interface
      interface
        subroutine gridInqYlongname(gridID,ylongname) bind(c,name='gridInqYlongname')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: ylongname
        end subroutine gridInqYlongname
      end interface
      interface
        subroutine gridDefYunits(gridID,yunits) bind(c,name='gridDefYunits')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: yunits
        end subroutine gridDefYunits
      end interface
      interface
        subroutine gridInqYunits(gridID,yunits) bind(c,name='gridInqYunits')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: yunits
        end subroutine gridInqYunits
      end interface
      interface
        subroutine gridInqXstdname(gridID,xstdname) bind(c,name='gridInqXstdname')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: xstdname
        end subroutine gridInqXstdname
      end interface
      interface
        subroutine gridInqYstdname(gridID,ystdname) bind(c,name='gridInqYstdname')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: ystdname
        end subroutine gridInqYstdname
      end interface
      interface
        subroutine gridDefPrec(gridID,prec) bind(c,name='gridDefPrec')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: prec
        end subroutine gridDefPrec
      end interface
      interface
        function gridInqPrec(gridID) bind(c,name='gridInqPrec')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqPrec
        end function gridInqPrec
      end interface
      interface
        function gridInqXval(gridID,index) bind(c,name='gridInqXval')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: index
          real(kind=c_double) :: gridInqXval
        end function gridInqXval
      end interface
      interface
        function gridInqYval(gridID,index) bind(c,name='gridInqYval')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: index
          real(kind=c_double) :: gridInqYval
        end function gridInqYval
      end interface
      interface
        function gridInqXinc(gridID) bind(c,name='gridInqXinc')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double) :: gridInqXinc
        end function gridInqXinc
      end interface
      interface
        function gridInqYinc(gridID) bind(c,name='gridInqYinc')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double) :: gridInqYinc
        end function gridInqYinc
      end interface
      interface
        function gridIsCircular(gridID) bind(c,name='gridIsCircular')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridIsCircular
        end function gridIsCircular
      end interface
      interface
        function gridIsRotated(gridID) bind(c,name='gridIsRotated')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridIsRotated
        end function gridIsRotated
      end interface
      interface
        subroutine gridDefXpole(gridID,xpole) bind(c,name='gridDefXpole')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), value :: xpole
        end subroutine gridDefXpole
      end interface
      interface
        function gridInqXpole(gridID) bind(c,name='gridInqXpole')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double) :: gridInqXpole
        end function gridInqXpole
      end interface
      interface
        subroutine gridDefYpole(gridID,ypole) bind(c,name='gridDefYpole')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), value :: ypole
        end subroutine gridDefYpole
      end interface
      interface
        function gridInqYpole(gridID) bind(c,name='gridInqYpole')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double) :: gridInqYpole
        end function gridInqYpole
      end interface
      interface
        subroutine gridDefAngle(gridID,angle) bind(c,name='gridDefAngle')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), value :: angle
        end subroutine gridDefAngle
      end interface
      interface
        function gridInqAngle(gridID) bind(c,name='gridInqAngle')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double) :: gridInqAngle
        end function gridInqAngle
      end interface
      interface
        function gridInqTrunc(gridID) bind(c,name='gridInqTrunc')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqTrunc
        end function gridInqTrunc
      end interface
      interface
        subroutine gridDefTrunc(gridID,trunc) bind(c,name='gridDefTrunc')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: trunc
        end subroutine gridDefTrunc
      end interface
      interface
        subroutine gridDefGMEnd(gridID,nd) bind(c,name='gridDefGMEnd')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: nd
        end subroutine gridDefGMEnd
      end interface
      interface
        function gridInqGMEnd(gridID) bind(c,name='gridInqGMEnd')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqGMEnd
        end function gridInqGMEnd
      end interface
      interface
        subroutine gridDefGMEni(gridID,ni) bind(c,name='gridDefGMEni')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: ni
        end subroutine gridDefGMEni
      end interface
      interface
        function gridInqGMEni(gridID) bind(c,name='gridInqGMEni')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqGMEni
        end function gridInqGMEni
      end interface
      interface
        subroutine gridDefGMEni2(gridID,ni2) bind(c,name='gridDefGMEni2')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: ni2
        end subroutine gridDefGMEni2
      end interface
      interface
        function gridInqGMEni2(gridID) bind(c,name='gridInqGMEni2')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqGMEni2
        end function gridInqGMEni2
      end interface
      interface
        subroutine gridDefGMEni3(gridID,ni3) bind(c,name='gridDefGMEni3')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: ni3
        end subroutine gridDefGMEni3
      end interface
      interface
        function gridInqGMEni3(gridID) bind(c,name='gridInqGMEni3')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqGMEni3
        end function gridInqGMEni3
      end interface
      interface
        subroutine gridDefNumber(gridID,number) bind(c,name='gridDefNumber')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: number
        end subroutine gridDefNumber
      end interface
      interface
        function gridInqNumber(gridID) bind(c,name='gridInqNumber')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqNumber
        end function gridInqNumber
      end interface
      interface
        subroutine gridDefPosition(gridID,position) bind(c,name='gridDefPosition')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: position
        end subroutine gridDefPosition
      end interface
      interface
        function gridInqPosition(gridID) bind(c,name='gridInqPosition')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqPosition
        end function gridInqPosition
      end interface
      interface
        subroutine gridDefReference(gridID,reference) bind(c,name='gridDefReference')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: reference
        end subroutine gridDefReference
      end interface
      interface
        function gridInqReference(gridID,reference) bind(c,name='gridInqReference')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(*) :: reference
          integer(kind=c_int) :: gridInqReference
        end function gridInqReference
      end interface
      interface
        subroutine gridDefUUID(gridID,uuid) bind(c,name='gridDefUUID')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(16) :: uuid
        end subroutine gridDefUUID
      end interface
      interface
        subroutine gridInqUUID(gridID,uuid) bind(c,name='gridInqUUID')
          import :: c_int,c_char
          integer(kind=c_int), value :: gridID
          character(kind=c_char), dimension(16) :: uuid
        end subroutine gridInqUUID
      end interface
      interface
        subroutine gridDefLCC(gridID,originLon,originLat,lonParY,lat1,lat2,xinc,yinc,projflag,scanflag) bind(c,name='gridDefLCC')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), value :: originLon
          real(kind=c_double), value :: originLat
          real(kind=c_double), value :: lonParY
          real(kind=c_double), value :: lat1
          real(kind=c_double), value :: lat2
          real(kind=c_double), value :: xinc
          real(kind=c_double), value :: yinc
          integer(kind=c_int), value :: projflag
          integer(kind=c_int), value :: scanflag
        end subroutine gridDefLCC
      end interface
      interface
        subroutine gridInqLCC(gridID,originLon,originLat,lonParY,lat1,lat2,xinc,yinc,projflag,scanflag) bind(c,name='gridInqLCC')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), intent(out) :: originLon
          real(kind=c_double), intent(out) :: originLat
          real(kind=c_double), intent(out) :: lonParY
          real(kind=c_double), intent(out) :: lat1
          real(kind=c_double), intent(out) :: lat2
          real(kind=c_double), intent(out) :: xinc
          real(kind=c_double), intent(out) :: yinc
          integer(kind=c_int), intent(out) :: projflag
          integer(kind=c_int), intent(out) :: scanflag
        end subroutine gridInqLCC
      end interface
      interface
        subroutine gridDefLcc2(gridID,earth_radius,lon_0,lat_0,lat_1,lat_2) bind(c,name='gridDefLcc2')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), value :: earth_radius
          real(kind=c_double), value :: lon_0
          real(kind=c_double), value :: lat_0
          real(kind=c_double), value :: lat_1
          real(kind=c_double), value :: lat_2
        end subroutine gridDefLcc2
      end interface
      interface
        subroutine gridInqLcc2(gridID,earth_radius,lon_0,lat_0,lat_1,lat_2) bind(c,name='gridInqLcc2')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), intent(out) :: earth_radius
          real(kind=c_double), intent(out) :: lon_0
          real(kind=c_double), intent(out) :: lat_0
          real(kind=c_double), intent(out) :: lat_1
          real(kind=c_double), intent(out) :: lat_2
        end subroutine gridInqLcc2
      end interface
      interface
        subroutine gridDefLaea(gridID,earth_radius,lon_0,lat_0) bind(c,name='gridDefLaea')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), value :: earth_radius
          real(kind=c_double), value :: lon_0
          real(kind=c_double), value :: lat_0
        end subroutine gridDefLaea
      end interface
      interface
        subroutine gridInqLaea(gridID,earth_radius,lon_0,lat_0) bind(c,name='gridInqLaea')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), intent(out) :: earth_radius
          real(kind=c_double), intent(out) :: lon_0
          real(kind=c_double), intent(out) :: lat_0
        end subroutine gridInqLaea
      end interface
      interface
        subroutine gridDefArea(gridID,area_vec) bind(c,name='gridDefArea')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), intent(in),dimension(*) :: area_vec
        end subroutine gridDefArea
      end interface
      interface
        subroutine gridInqArea(gridID,area_vec) bind(c,name='gridInqArea')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), intent(out),dimension(*) :: area_vec
        end subroutine gridInqArea
      end interface
      interface
        function gridHasArea(gridID) bind(c,name='gridHasArea')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridHasArea
        end function gridHasArea
      end interface
      interface
        subroutine gridDefNvertex(gridID,nvertex) bind(c,name='gridDefNvertex')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: nvertex
        end subroutine gridDefNvertex
      end interface
      interface
        function gridInqNvertex(gridID) bind(c,name='gridInqNvertex')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqNvertex
        end function gridInqNvertex
      end interface
      interface
        subroutine gridDefXbounds(gridID,xbounds_vec) bind(c,name='gridDefXbounds')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), intent(in),dimension(*) :: xbounds_vec
        end subroutine gridDefXbounds
      end interface
      interface
        function gridInqXbounds(gridID,xbounds_vec) bind(c,name='gridInqXbounds')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), intent(out),dimension(*) :: xbounds_vec
          integer(kind=c_int) :: gridInqXbounds
        end function gridInqXbounds
      end interface
      interface
        subroutine gridDefYbounds(gridID,ybounds_vec) bind(c,name='gridDefYbounds')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), intent(in),dimension(*) :: ybounds_vec
        end subroutine gridDefYbounds
      end interface
      interface
        function gridInqYbounds(gridID,ybounds_vec) bind(c,name='gridInqYbounds')
          import :: c_int,c_double
          integer(kind=c_int), value :: gridID
          real(kind=c_double), intent(out),dimension(*) :: ybounds_vec
          integer(kind=c_int) :: gridInqYbounds
        end function gridInqYbounds
      end interface
      interface
        subroutine gridDefRowlon(gridID,nrowlon,rowlon_vec) bind(c,name='gridDefRowlon')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: nrowlon
          integer(kind=c_int), intent(in),dimension(*) :: rowlon_vec
        end subroutine gridDefRowlon
      end interface
      interface
        subroutine gridInqRowlon(gridID,rowlon_vec) bind(c,name='gridInqRowlon')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), intent(out),dimension(*) :: rowlon_vec
        end subroutine gridInqRowlon
      end interface
      interface
        subroutine gridChangeType(gridID,gridtype) bind(c,name='gridChangeType')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: gridtype
        end subroutine gridChangeType
      end interface
      interface
        subroutine gridDefComplexPacking(gridID,lpack) bind(c,name='gridDefComplexPacking')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: lpack
        end subroutine gridDefComplexPacking
      end interface
      interface
        function gridInqComplexPacking(gridID) bind(c,name='gridInqComplexPacking')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: gridInqComplexPacking
        end function gridInqComplexPacking
      end interface
      interface
        subroutine zaxisName(zaxistype,zaxisnamev) bind(c,name='zaxisName')
          import :: c_int,c_char
          integer(kind=c_int), value :: zaxistype
          character(kind=c_char), dimension(*) :: zaxisnamev
        end subroutine zaxisName
      end interface
      interface
        function zaxisCreate(zaxistype,size) bind(c,name='zaxisCreate')
          import :: c_int
          integer(kind=c_int), value :: zaxistype
          integer(kind=c_int), value :: size
          integer(kind=c_int) :: zaxisCreate
        end function zaxisCreate
      end interface
      interface
        subroutine zaxisDestroy(zaxisID) bind(c,name='zaxisDestroy')
          import :: c_int
          integer(kind=c_int), value :: zaxisID
        end subroutine zaxisDestroy
      end interface
      interface
        function zaxisInqType(zaxisID) bind(c,name='zaxisInqType')
          import :: c_int
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int) :: zaxisInqType
        end function zaxisInqType
      end interface
      interface
        function zaxisInqSize(zaxisID) bind(c,name='zaxisInqSize')
          import :: c_int
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int) :: zaxisInqSize
        end function zaxisInqSize
      end interface
      interface
        function zaxisDuplicate(zaxisID) bind(c,name='zaxisDuplicate')
          import :: c_int
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int) :: zaxisDuplicate
        end function zaxisDuplicate
      end interface
      interface
        subroutine zaxisResize(zaxisID,size) bind(c,name='zaxisResize')
          import :: c_int
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int), value :: size
        end subroutine zaxisResize
      end interface
      interface
        subroutine zaxisPrint(zaxisID) bind(c,name='zaxisPrint')
          import :: c_int
          integer(kind=c_int), value :: zaxisID
        end subroutine zaxisPrint
      end interface
      interface
        subroutine zaxisDefLevels(zaxisID,levels_vec) bind(c,name='zaxisDefLevels')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          real(kind=c_double), intent(in),dimension(*) :: levels_vec
        end subroutine zaxisDefLevels
      end interface
      interface
        subroutine zaxisInqLevels(zaxisID,levels_vec) bind(c,name='zaxisInqLevels')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          real(kind=c_double), intent(out),dimension(*) :: levels_vec
        end subroutine zaxisInqLevels
      end interface
      interface
        subroutine zaxisDefLevel(zaxisID,levelID,levels) bind(c,name='zaxisDefLevel')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int), value :: levelID
          real(kind=c_double), value :: levels
        end subroutine zaxisDefLevel
      end interface
      interface
        function zaxisInqLevel(zaxisID,levelID) bind(c,name='zaxisInqLevel')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int), value :: levelID
          real(kind=c_double) :: zaxisInqLevel
        end function zaxisInqLevel
      end interface
      interface
        subroutine zaxisDefNlevRef(gridID,nhlev) bind(c,name='zaxisDefNlevRef')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: nhlev
        end subroutine zaxisDefNlevRef
      end interface
      interface
        function zaxisInqNlevRef(gridID) bind(c,name='zaxisInqNlevRef')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: zaxisInqNlevRef
        end function zaxisInqNlevRef
      end interface
      interface
        subroutine zaxisDefNumber(gridID,number) bind(c,name='zaxisDefNumber')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int), value :: number
        end subroutine zaxisDefNumber
      end interface
      interface
        function zaxisInqNumber(gridID) bind(c,name='zaxisInqNumber')
          import :: c_int
          integer(kind=c_int), value :: gridID
          integer(kind=c_int) :: zaxisInqNumber
        end function zaxisInqNumber
      end interface
      interface
        subroutine zaxisDefUUID(zaxisID,uuid) bind(c,name='zaxisDefUUID')
          import :: c_int,c_char
          integer(kind=c_int), value :: zaxisID
          character(kind=c_char), dimension(16) :: uuid
        end subroutine zaxisDefUUID
      end interface
      interface
        subroutine zaxisInqUUID(zaxisID,uuid) bind(c,name='zaxisInqUUID')
          import :: c_int,c_char
          integer(kind=c_int), value :: zaxisID
          character(kind=c_char), dimension(16) :: uuid
        end subroutine zaxisInqUUID
      end interface
      interface
        subroutine zaxisDefName(zaxisID,name) bind(c,name='zaxisDefName')
          import :: c_int,c_char
          integer(kind=c_int), value :: zaxisID
          character(kind=c_char), dimension(*) :: name
        end subroutine zaxisDefName
      end interface
      interface
        subroutine zaxisInqName(zaxisID,name) bind(c,name='zaxisInqName')
          import :: c_int,c_char
          integer(kind=c_int), value :: zaxisID
          character(kind=c_char), dimension(*) :: name
        end subroutine zaxisInqName
      end interface
      interface
        subroutine zaxisDefLongname(zaxisID,longname) bind(c,name='zaxisDefLongname')
          import :: c_int,c_char
          integer(kind=c_int), value :: zaxisID
          character(kind=c_char), dimension(*) :: longname
        end subroutine zaxisDefLongname
      end interface
      interface
        subroutine zaxisInqLongname(zaxisID,longname) bind(c,name='zaxisInqLongname')
          import :: c_int,c_char
          integer(kind=c_int), value :: zaxisID
          character(kind=c_char), dimension(*) :: longname
        end subroutine zaxisInqLongname
      end interface
      interface
        subroutine zaxisDefUnits(zaxisID,units) bind(c,name='zaxisDefUnits')
          import :: c_int,c_char
          integer(kind=c_int), value :: zaxisID
          character(kind=c_char), dimension(*) :: units
        end subroutine zaxisDefUnits
      end interface
      interface
        subroutine zaxisInqUnits(zaxisID,units) bind(c,name='zaxisInqUnits')
          import :: c_int,c_char
          integer(kind=c_int), value :: zaxisID
          character(kind=c_char), dimension(*) :: units
        end subroutine zaxisInqUnits
      end interface
      interface
        subroutine zaxisInqStdname(zaxisID,stdname) bind(c,name='zaxisInqStdname')
          import :: c_int,c_char
          integer(kind=c_int), value :: zaxisID
          character(kind=c_char), dimension(*) :: stdname
        end subroutine zaxisInqStdname
      end interface
      interface
        subroutine zaxisDefPrec(zaxisID,prec) bind(c,name='zaxisDefPrec')
          import :: c_int
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int), value :: prec
        end subroutine zaxisDefPrec
      end interface
      interface
        function zaxisInqPrec(zaxisID) bind(c,name='zaxisInqPrec')
          import :: c_int
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int) :: zaxisInqPrec
        end function zaxisInqPrec
      end interface
      interface
        subroutine zaxisDefPositive(zaxisID,positive) bind(c,name='zaxisDefPositive')
          import :: c_int
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int), value :: positive
        end subroutine zaxisDefPositive
      end interface
      interface
        function zaxisInqPositive(zaxisID) bind(c,name='zaxisInqPositive')
          import :: c_int
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int) :: zaxisInqPositive
        end function zaxisInqPositive
      end interface
      interface
        subroutine zaxisDefLtype(zaxisID,ltype) bind(c,name='zaxisDefLtype')
          import :: c_int
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int), value :: ltype
        end subroutine zaxisDefLtype
      end interface
      interface
        function zaxisInqLtype(zaxisID) bind(c,name='zaxisInqLtype')
          import :: c_int
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int) :: zaxisInqLtype
        end function zaxisInqLtype
      end interface
      interface
        function zaxisInqLevelsPtr(zaxisID) bind(c,name='zaxisInqLevelsPtr')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          real(kind=c_double) :: zaxisInqLevelsPtr
        end function zaxisInqLevelsPtr
      end interface
      interface
        subroutine zaxisDefVct(zaxisID,size,vct_vec) bind(c,name='zaxisDefVct')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int), value :: size
          real(kind=c_double), intent(in),dimension(*) :: vct_vec
        end subroutine zaxisDefVct
      end interface
      interface
        subroutine zaxisInqVct(zaxisID,vct_vec) bind(c,name='zaxisInqVct')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          real(kind=c_double), intent(out),dimension(*) :: vct_vec
        end subroutine zaxisInqVct
      end interface
      interface
        function zaxisInqVctSize(zaxisID) bind(c,name='zaxisInqVctSize')
          import :: c_int
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int) :: zaxisInqVctSize
        end function zaxisInqVctSize
      end interface
      interface
        function zaxisInqVctPtr(zaxisID) bind(c,name='zaxisInqVctPtr')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          real(kind=c_double) :: zaxisInqVctPtr
        end function zaxisInqVctPtr
      end interface
      interface
        subroutine zaxisDefLbounds(zaxisID,lbounds_vec) bind(c,name='zaxisDefLbounds')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          real(kind=c_double), intent(in),dimension(*) :: lbounds_vec
        end subroutine zaxisDefLbounds
      end interface
      interface
        function zaxisInqLbounds(zaxisID,lbounds_vec) bind(c,name='zaxisInqLbounds')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          real(kind=c_double), intent(out),dimension(*) :: lbounds_vec
          integer(kind=c_int) :: zaxisInqLbounds
        end function zaxisInqLbounds
      end interface
      interface
        function zaxisInqLbound(zaxisID,index) bind(c,name='zaxisInqLbound')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int), value :: index
          real(kind=c_double) :: zaxisInqLbound
        end function zaxisInqLbound
      end interface
      interface
        subroutine zaxisDefUbounds(zaxisID,ubounds_vec) bind(c,name='zaxisDefUbounds')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          real(kind=c_double), intent(in),dimension(*) :: ubounds_vec
        end subroutine zaxisDefUbounds
      end interface
      interface
        function zaxisInqUbounds(zaxisID,ubounds_vec) bind(c,name='zaxisInqUbounds')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          real(kind=c_double), intent(out),dimension(*) :: ubounds_vec
          integer(kind=c_int) :: zaxisInqUbounds
        end function zaxisInqUbounds
      end interface
      interface
        function zaxisInqUbound(zaxisID,index) bind(c,name='zaxisInqUbound')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int), value :: index
          real(kind=c_double) :: zaxisInqUbound
        end function zaxisInqUbound
      end interface
      interface
        subroutine zaxisDefWeights(zaxisID,weights_vec) bind(c,name='zaxisDefWeights')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          real(kind=c_double), intent(in),dimension(*) :: weights_vec
        end subroutine zaxisDefWeights
      end interface
      interface
        function zaxisInqWeights(zaxisID,weights_vec) bind(c,name='zaxisInqWeights')
          import :: c_int,c_double
          integer(kind=c_int), value :: zaxisID
          real(kind=c_double), intent(out),dimension(*) :: weights_vec
          integer(kind=c_int) :: zaxisInqWeights
        end function zaxisInqWeights
      end interface
      interface
        subroutine zaxisChangeType(zaxisID,zaxistype) bind(c,name='zaxisChangeType')
          import :: c_int
          integer(kind=c_int), value :: zaxisID
          integer(kind=c_int), value :: zaxistype
        end subroutine zaxisChangeType
      end interface
      interface
        function taxisCreate(timetype) bind(c,name='taxisCreate')
          import :: c_int
          integer(kind=c_int), value :: timetype
          integer(kind=c_int) :: taxisCreate
        end function taxisCreate
      end interface
      interface
        subroutine taxisDestroy(taxisID) bind(c,name='taxisDestroy')
          import :: c_int
          integer(kind=c_int), value :: taxisID
        end subroutine taxisDestroy
      end interface
      interface
        function taxisDuplicate(taxisID) bind(c,name='taxisDuplicate')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int) :: taxisDuplicate
        end function taxisDuplicate
      end interface
      interface
        subroutine taxisCopyTimestep(taxisIDdes,taxisIDsrc) bind(c,name='taxisCopyTimestep')
          import :: c_int
          integer(kind=c_int), value :: taxisIDdes
          integer(kind=c_int), value :: taxisIDsrc
        end subroutine taxisCopyTimestep
      end interface
      interface
        subroutine taxisDefType(taxisID,type) bind(c,name='taxisDefType')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), value :: type
        end subroutine taxisDefType
      end interface
      interface
        subroutine taxisDefVdate(taxisID,date) bind(c,name='taxisDefVdate')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), value :: date
        end subroutine taxisDefVdate
      end interface
      interface
        subroutine taxisDefVtime(taxisID,time) bind(c,name='taxisDefVtime')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), value :: time
        end subroutine taxisDefVtime
      end interface
      interface
        function taxisInqVdate(taxisID) bind(c,name='taxisInqVdate')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int) :: taxisInqVdate
        end function taxisInqVdate
      end interface
      interface
        function taxisInqVtime(taxisID) bind(c,name='taxisInqVtime')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int) :: taxisInqVtime
        end function taxisInqVtime
      end interface
      interface
        subroutine taxisDefRdate(taxisID,date) bind(c,name='taxisDefRdate')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), value :: date
        end subroutine taxisDefRdate
      end interface
      interface
        subroutine taxisDefRtime(taxisID,time) bind(c,name='taxisDefRtime')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), value :: time
        end subroutine taxisDefRtime
      end interface
      interface
        function taxisInqRdate(taxisID) bind(c,name='taxisInqRdate')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int) :: taxisInqRdate
        end function taxisInqRdate
      end interface
      interface
        function taxisInqRtime(taxisID) bind(c,name='taxisInqRtime')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int) :: taxisInqRtime
        end function taxisInqRtime
      end interface
      interface
        subroutine taxisDefFdate(taxisID,date) bind(c,name='taxisDefFdate')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), value :: date
        end subroutine taxisDefFdate
      end interface
      interface
        subroutine taxisDefFtime(taxisID,time) bind(c,name='taxisDefFtime')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), value :: time
        end subroutine taxisDefFtime
      end interface
      interface
        function taxisInqFdate(taxisID) bind(c,name='taxisInqFdate')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int) :: taxisInqFdate
        end function taxisInqFdate
      end interface
      interface
        function taxisInqFtime(taxisID) bind(c,name='taxisInqFtime')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int) :: taxisInqFtime
        end function taxisInqFtime
      end interface
      interface
        function taxisHasBounds(taxisID) bind(c,name='taxisHasBounds')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int) :: taxisHasBounds
        end function taxisHasBounds
      end interface
      interface
        subroutine taxisDeleteBounds(taxisID) bind(c,name='taxisDeleteBounds')
          import :: c_int
          integer(kind=c_int), value :: taxisID
        end subroutine taxisDeleteBounds
      end interface
      interface
        subroutine taxisDefVdateBounds(taxisID,vdate_lb,vdate_ub) bind(c,name='taxisDefVdateBounds')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), value :: vdate_lb
          integer(kind=c_int), value :: vdate_ub
        end subroutine taxisDefVdateBounds
      end interface
      interface
        subroutine taxisDefVtimeBounds(taxisID,vtime_lb,vtime_ub) bind(c,name='taxisDefVtimeBounds')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), value :: vtime_lb
          integer(kind=c_int), value :: vtime_ub
        end subroutine taxisDefVtimeBounds
      end interface
      interface
        subroutine taxisInqVdateBounds(taxisID,vdate_lb,vdate_ub) bind(c,name='taxisInqVdateBounds')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), intent(out) :: vdate_lb
          integer(kind=c_int), intent(out) :: vdate_ub
        end subroutine taxisInqVdateBounds
      end interface
      interface
        subroutine taxisInqVtimeBounds(taxisID,vtime_lb,vtime_ub) bind(c,name='taxisInqVtimeBounds')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), intent(out) :: vtime_lb
          integer(kind=c_int), intent(out) :: vtime_ub
        end subroutine taxisInqVtimeBounds
      end interface
      interface
        subroutine taxisDefCalendar(taxisID,calendar) bind(c,name='taxisDefCalendar')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), value :: calendar
        end subroutine taxisDefCalendar
      end interface
      interface
        function taxisInqCalendar(taxisID) bind(c,name='taxisInqCalendar')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int) :: taxisInqCalendar
        end function taxisInqCalendar
      end interface
      interface
        subroutine taxisDefTunit(taxisID,tunit) bind(c,name='taxisDefTunit')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), value :: tunit
        end subroutine taxisDefTunit
      end interface
      interface
        function taxisInqTunit(taxisID) bind(c,name='taxisInqTunit')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int) :: taxisInqTunit
        end function taxisInqTunit
      end interface
      interface
        subroutine taxisDefForecastTunit(taxisID,tunit) bind(c,name='taxisDefForecastTunit')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), value :: tunit
        end subroutine taxisDefForecastTunit
      end interface
      interface
        function taxisInqForecastTunit(taxisID) bind(c,name='taxisInqForecastTunit')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int) :: taxisInqForecastTunit
        end function taxisInqForecastTunit
      end interface
      interface
        subroutine taxisDefForecastPeriod(taxisID,fc_period) bind(c,name='taxisDefForecastPeriod')
          import :: c_int,c_double
          integer(kind=c_int), value :: taxisID
          real(kind=c_double), value :: fc_period
        end subroutine taxisDefForecastPeriod
      end interface
      interface
        function taxisInqForecastPeriod(taxisID) bind(c,name='taxisInqForecastPeriod')
          import :: c_int,c_double
          integer(kind=c_int), value :: taxisID
          real(kind=c_double) :: taxisInqForecastPeriod
        end function taxisInqForecastPeriod
      end interface
      interface
        subroutine taxisDefNumavg(taxisID,numavg) bind(c,name='taxisDefNumavg')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int), value :: numavg
        end subroutine taxisDefNumavg
      end interface
      interface
        function taxisInqType(taxisID) bind(c,name='taxisInqType')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int) :: taxisInqType
        end function taxisInqType
      end interface
      interface
        function taxisInqNumavg(taxisID) bind(c,name='taxisInqNumavg')
          import :: c_int
          integer(kind=c_int), value :: taxisID
          integer(kind=c_int) :: taxisInqNumavg
        end function taxisInqNumavg
      end interface
      interface
        function institutDef(center,subcenter,name,longname) bind(c,name='institutDef')
          import :: c_int,c_char
          integer(kind=c_int), value :: center
          integer(kind=c_int), value :: subcenter
          character(kind=c_char), dimension(*) :: name
          character(kind=c_char), dimension(*) :: longname
          integer(kind=c_int) :: institutDef
        end function institutDef
      end interface
      interface
        function institutInq(center,subcenter,name,longname) bind(c,name='institutInq')
          import :: c_int,c_char
          integer(kind=c_int), value :: center
          integer(kind=c_int), value :: subcenter
          character(kind=c_char), dimension(*) :: name
          character(kind=c_char), dimension(*) :: longname
          integer(kind=c_int) :: institutInq
        end function institutInq
      end interface
      interface
        function institutInqNumber() bind(c,name='institutInqNumber')
          import :: c_int
          integer(kind=c_int) :: institutInqNumber
        end function institutInqNumber
      end interface
      interface
        function institutInqCenter(instID) bind(c,name='institutInqCenter')
          import :: c_int
          integer(kind=c_int), value :: instID
          integer(kind=c_int) :: institutInqCenter
        end function institutInqCenter
      end interface
      interface
        function institutInqSubcenter(instID) bind(c,name='institutInqSubcenter')
          import :: c_int
          integer(kind=c_int), value :: instID
          integer(kind=c_int) :: institutInqSubcenter
        end function institutInqSubcenter
      end interface
      interface
        function modelDef(instID,modelgribID,name) bind(c,name='modelDef')
          import :: c_int,c_char
          integer(kind=c_int), value :: instID
          integer(kind=c_int), value :: modelgribID
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int) :: modelDef
        end function modelDef
      end interface
      interface
        function modelInq(instID,modelgribID,name) bind(c,name='modelInq')
          import :: c_int,c_char
          integer(kind=c_int), value :: instID
          integer(kind=c_int), value :: modelgribID
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int) :: modelInq
        end function modelInq
      end interface
      interface
        function modelInqInstitut(modelID) bind(c,name='modelInqInstitut')
          import :: c_int
          integer(kind=c_int), value :: modelID
          integer(kind=c_int) :: modelInqInstitut
        end function modelInqInstitut
      end interface
      interface
        function modelInqGribID(modelID) bind(c,name='modelInqGribID')
          import :: c_int
          integer(kind=c_int), value :: modelID
          integer(kind=c_int) :: modelInqGribID
        end function modelInqGribID
      end interface
      interface
        subroutine tableWriteC(filename,tableID) bind(c,name='tableWriteC')
          import :: c_char,c_int
          character(kind=c_char), dimension(*) :: filename
          integer(kind=c_int), value :: tableID
        end subroutine tableWriteC
      end interface
      interface
        subroutine tableWrite(filename,tableID) bind(c,name='tableWrite')
          import :: c_char,c_int
          character(kind=c_char), dimension(*) :: filename
          integer(kind=c_int), value :: tableID
        end subroutine tableWrite
      end interface
      interface
        function tableRead(tablefile) bind(c,name='tableRead')
          import :: c_char,c_int
          character(kind=c_char), dimension(*) :: tablefile
          integer(kind=c_int) :: tableRead
        end function tableRead
      end interface
      interface
        function tableDef(modelID,tablenum,tablename) bind(c,name='tableDef')
          import :: c_int,c_char
          integer(kind=c_int), value :: modelID
          integer(kind=c_int), value :: tablenum
          character(kind=c_char), dimension(*) :: tablename
          integer(kind=c_int) :: tableDef
        end function tableDef
      end interface
      interface
        subroutine tableDefEntry(tableID,code,name,longname,units) bind(c,name='tableDefEntry')
          import :: c_int,c_char
          integer(kind=c_int), value :: tableID
          integer(kind=c_int), value :: code
          character(kind=c_char), dimension(*) :: name
          character(kind=c_char), dimension(*) :: longname
          character(kind=c_char), dimension(*) :: units
        end subroutine tableDefEntry
      end interface
      interface
        function tableInq(modelID,tablenum,tablename) bind(c,name='tableInq')
          import :: c_int,c_char
          integer(kind=c_int), value :: modelID
          integer(kind=c_int), value :: tablenum
          character(kind=c_char), dimension(*) :: tablename
          integer(kind=c_int) :: tableInq
        end function tableInq
      end interface
      interface
        function tableInqNumber() bind(c,name='tableInqNumber')
          import :: c_int
          integer(kind=c_int) :: tableInqNumber
        end function tableInqNumber
      end interface
      interface
        function tableInqNum(tableID) bind(c,name='tableInqNum')
          import :: c_int
          integer(kind=c_int), value :: tableID
          integer(kind=c_int) :: tableInqNum
        end function tableInqNum
      end interface
      interface
        function tableInqModel(tableID) bind(c,name='tableInqModel')
          import :: c_int
          integer(kind=c_int), value :: tableID
          integer(kind=c_int) :: tableInqModel
        end function tableInqModel
      end interface
      interface
        subroutine tableInqPar(tableID,code,name,longname,units) bind(c,name='tableInqPar')
          import :: c_int,c_char
          integer(kind=c_int), value :: tableID
          integer(kind=c_int), value :: code
          character(kind=c_char), dimension(*) :: name
          character(kind=c_char), dimension(*) :: longname
          character(kind=c_char), dimension(*) :: units
        end subroutine tableInqPar
      end interface
      interface
        function tableInqParCode(tableID,name,code) bind(c,name='tableInqParCode')
          import :: c_int,c_char
          integer(kind=c_int), value :: tableID
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int), intent(out) :: code
          integer(kind=c_int) :: tableInqParCode
        end function tableInqParCode
      end interface
      interface
        function tableInqParName(tableID,code,name) bind(c,name='tableInqParName')
          import :: c_int,c_char
          integer(kind=c_int), value :: tableID
          integer(kind=c_int), value :: code
          character(kind=c_char), dimension(*) :: name
          integer(kind=c_int) :: tableInqParName
        end function tableInqParName
      end interface
      interface
        function tableInqParLongname(tableID,code,longname) bind(c,name='tableInqParLongname')
          import :: c_int,c_char
          integer(kind=c_int), value :: tableID
          integer(kind=c_int), value :: code
          character(kind=c_char), dimension(*) :: longname
          integer(kind=c_int) :: tableInqParLongname
        end function tableInqParLongname
      end interface
      interface
        function tableInqParUnits(tableID,code,units) bind(c,name='tableInqParUnits')
          import :: c_int,c_char
          integer(kind=c_int), value :: tableID
          integer(kind=c_int), value :: code
          character(kind=c_char), dimension(*) :: units
          integer(kind=c_int) :: tableInqParUnits
        end function tableInqParUnits
      end interface
      interface
        subroutine streamDefHistory(streamID,size,history) bind(c,name='streamDefHistory')
          import :: c_int,c_char
          integer(kind=c_int), value :: streamID
          integer(kind=c_int), value :: size
          character(kind=c_char), dimension(*) :: history
        end subroutine streamDefHistory
      end interface
      interface
        function streamInqHistorySize(streamID) bind(c,name='streamInqHistorySize')
          import :: c_int
          integer(kind=c_int), value :: streamID
          integer(kind=c_int) :: streamInqHistorySize
        end function streamInqHistorySize
      end interface
      interface
        subroutine streamInqHistoryString(streamID,history) bind(c,name='streamInqHistoryString')
          import :: c_int,c_char
          integer(kind=c_int), value :: streamID
          character(kind=c_char), dimension(*) :: history
        end subroutine streamInqHistoryString
      end interface
      interface
        subroutine gribapiLibraryVersion(major_version,minor_version,revision_version) bind(c,name='gribapiLibraryVersion')
          import :: c_int
          integer(kind=c_int), intent(out) :: major_version
          integer(kind=c_int), intent(out) :: minor_version
          integer(kind=c_int), intent(out) :: revision_version
        end subroutine gribapiLibraryVersion
      end interface

      public :: strlen
      public :: getchar
      public :: getchar_unlocked
      public :: cdiReset
      public :: cdiStringError
      public :: cdiDebug
      public :: cdiLibraryVersion
      public :: cdiPrintVersion
      public :: cdiHaveFiletype
      public :: cdiDefMissval
      public :: cdiInqMissval
      public :: cdiDefGlobal
      public :: namespaceNew
      public :: namespaceSetActive
      public :: namespaceDelete
      public :: cdiParamToString
      public :: cdiDecodeParam
      public :: cdiEncodeParam
      public :: cdiDecodeDate
      public :: cdiEncodeDate
      public :: cdiDecodeTime
      public :: cdiEncodeTime
      public :: cdiGetFiletype
      public :: streamOpenRead
      public :: streamOpenWrite
      public :: streamOpenAppend
      public :: streamClose
      public :: streamSync
      public :: streamDefVlist
      public :: streamInqVlist
      public :: streamInqVlistIDorig
      public :: streamInqFiletype
      public :: streamDefByteorder
      public :: streamInqByteorder
      public :: streamDefCompType
      public :: streamInqCompType
      public :: streamDefCompLevel
      public :: streamInqCompLevel
      public :: streamDefTimestep
      public :: streamInqTimestep
      public :: streamInqCurTimestepID
      public :: streamFilename
      public :: streamFilesuffix
      public :: streamInqNvars
      public :: streamWriteVar
      public :: streamWriteVarF
      public :: streamReadVar
      public :: streamReadVarF
      public :: streamWriteVarSlice
      public :: streamWriteVarSliceF
      public :: streamReadVarSlice
      public :: streamReadVarSliceF
      public :: streamDefRecord
      public :: streamInqRecord
      public :: streamWriteRecord
      public :: streamWriteRecordF
      public :: streamReadRecord
      public :: streamCopyRecord
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
      public :: vlistInqVarID
      public :: vlistDefVarTsteptype
      public :: vlistInqVarTsteptype
      public :: vlistDefVarCompType
      public :: vlistInqVarCompType
      public :: vlistDefVarCompLevel
      public :: vlistInqVarCompLevel
      public :: vlistDefVarParam
      public :: vlistInqVarParam
      public :: vlistDefVarCode
      public :: vlistInqVarCode
      public :: vlistDefVarDatatype
      public :: vlistInqVarDatatype
      public :: vlistDefVarChunkType
      public :: vlistInqVarChunkType
      public :: vlistDefVarXYZ
      public :: vlistInqVarXYZ
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
      public :: vlistDefVarExtra
      public :: vlistInqVarExtra
      public :: vlistDefVarScalefactor
      public :: vlistInqVarScalefactor
      public :: vlistDefVarAddoffset
      public :: vlistInqVarAddoffset
      public :: vlistDefVarTimave
      public :: vlistInqVarTimave
      public :: vlistDefVarTimaccu
      public :: vlistInqVarTimaccu
      public :: vlistDefVarTypeOfGeneratingProcess
      public :: vlistInqVarTypeOfGeneratingProcess
      public :: vlistDefVarProductDefinitionTemplate
      public :: vlistInqVarProductDefinitionTemplate
      public :: vlistInqVarSize
      public :: vlistDefIndex
      public :: vlistInqIndex
      public :: vlistDefFlag
      public :: vlistInqFlag
      public :: vlistFindVar
      public :: vlistFindLevel
      public :: vlistMergedVar
      public :: vlistMergedLevel
      public :: vlistDefVarEnsemble
      public :: vlistInqVarEnsemble
      public :: cdiClearAdditionalKeys
      public :: cdiDefAdditionalKey
      public :: vlistDefVarIntKey
      public :: vlistDefVarDblKey
      public :: vlistHasVarKey
      public :: vlistInqVarDblKey
      public :: vlistInqVarIntKey
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
      public :: gridDefMaskGME
      public :: gridInqMaskGME
      public :: gridDefMask
      public :: gridInqMask
      public :: gridPrint
      public :: gridCreate
      public :: gridDestroy
      public :: gridDuplicate
      public :: gridInqType
      public :: gridInqSize
      public :: gridDefXsize
      public :: gridInqXsize
      public :: gridDefYsize
      public :: gridInqYsize
      public :: gridDefNP
      public :: gridInqNP
      public :: gridDefXvals
      public :: gridInqXvals
      public :: gridDefYvals
      public :: gridInqYvals
      public :: gridDefXname
      public :: gridInqXname
      public :: gridDefXlongname
      public :: gridInqXlongname
      public :: gridDefXunits
      public :: gridInqXunits
      public :: gridDefYname
      public :: gridInqYname
      public :: gridDefYlongname
      public :: gridInqYlongname
      public :: gridDefYunits
      public :: gridInqYunits
      public :: gridInqXstdname
      public :: gridInqYstdname
      public :: gridDefPrec
      public :: gridInqPrec
      public :: gridInqXval
      public :: gridInqYval
      public :: gridInqXinc
      public :: gridInqYinc
      public :: gridIsCircular
      public :: gridIsRotated
      public :: gridDefXpole
      public :: gridInqXpole
      public :: gridDefYpole
      public :: gridInqYpole
      public :: gridDefAngle
      public :: gridInqAngle
      public :: gridInqTrunc
      public :: gridDefTrunc
      public :: gridDefGMEnd
      public :: gridInqGMEnd
      public :: gridDefGMEni
      public :: gridInqGMEni
      public :: gridDefGMEni2
      public :: gridInqGMEni2
      public :: gridDefGMEni3
      public :: gridInqGMEni3
      public :: gridDefNumber
      public :: gridInqNumber
      public :: gridDefPosition
      public :: gridInqPosition
      public :: gridDefReference
      public :: gridInqReference
      public :: gridDefUUID
      public :: gridInqUUID
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
      public :: zaxisDefLevels
      public :: zaxisInqLevels
      public :: zaxisDefLevel
      public :: zaxisInqLevel
      public :: zaxisDefNlevRef
      public :: zaxisInqNlevRef
      public :: zaxisDefNumber
      public :: zaxisInqNumber
      public :: zaxisDefUUID
      public :: zaxisInqUUID
      public :: zaxisDefName
      public :: zaxisInqName
      public :: zaxisDefLongname
      public :: zaxisInqLongname
      public :: zaxisDefUnits
      public :: zaxisInqUnits
      public :: zaxisInqStdname
      public :: zaxisDefPrec
      public :: zaxisInqPrec
      public :: zaxisDefPositive
      public :: zaxisInqPositive
      public :: zaxisDefLtype
      public :: zaxisInqLtype
      public :: zaxisInqLevelsPtr
      public :: zaxisDefVct
      public :: zaxisInqVct
      public :: zaxisInqVctSize
      public :: zaxisInqVctPtr
      public :: zaxisDefLbounds
      public :: zaxisInqLbounds
      public :: zaxisInqLbound
      public :: zaxisDefUbounds
      public :: zaxisInqUbounds
      public :: zaxisInqUbound
      public :: zaxisDefWeights
      public :: zaxisInqWeights
      public :: zaxisChangeType
      public :: taxisCreate
      public :: taxisDestroy
      public :: taxisDuplicate
      public :: taxisCopyTimestep
      public :: taxisDefType
      public :: taxisDefVdate
      public :: taxisDefVtime
      public :: taxisInqVdate
      public :: taxisInqVtime
      public :: taxisDefRdate
      public :: taxisDefRtime
      public :: taxisInqRdate
      public :: taxisInqRtime
      public :: taxisDefFdate
      public :: taxisDefFtime
      public :: taxisInqFdate
      public :: taxisInqFtime
      public :: taxisHasBounds
      public :: taxisDeleteBounds
      public :: taxisDefVdateBounds
      public :: taxisDefVtimeBounds
      public :: taxisInqVdateBounds
      public :: taxisInqVtimeBounds
      public :: taxisDefCalendar
      public :: taxisInqCalendar
      public :: taxisDefTunit
      public :: taxisInqTunit
      public :: taxisDefForecastTunit
      public :: taxisInqForecastTunit
      public :: taxisDefForecastPeriod
      public :: taxisInqForecastPeriod
      public :: taxisDefNumavg
      public :: taxisInqType
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
      public :: gribapiLibraryVersion
      public :: ctrim
      public :: c_len

      public :: CDI_MAX_NAME
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
      public :: FILETYPE_UNDEF
      public :: FILETYPE_GRB
      public :: FILETYPE_GRB2
      public :: FILETYPE_NC
      public :: FILETYPE_NC2
      public :: FILETYPE_NC4
      public :: FILETYPE_NC4C
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
      public :: DATATYPE_UCHAR
      public :: CHUNK_AUTO
      public :: CHUNK_GRID
      public :: CHUNK_LINES
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
      public :: GRID_PROJECTION
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
      public :: ZAXIS_TOA
      public :: ZAXIS_SEA_BOTTOM
      public :: ZAXIS_ATMOSPHERE
      public :: ZAXIS_CLOUD_BASE
      public :: ZAXIS_CLOUD_TOP
      public :: ZAXIS_ISOTHERM_ZERO
      public :: ZAXIS_SNOW
      public :: ZAXIS_LAKE_BOTTOM
      public :: ZAXIS_SEDIMENT_BOTTOM
      public :: ZAXIS_SEDIMENT_BOTTOM_TA
      public :: ZAXIS_SEDIMENT_BOTTOM_TW
      public :: ZAXIS_MIX_LAYER
      public :: ZAXIS_REFERENCE
      public :: TIME_CONSTANT
      public :: TIME_VARIABLE
      public :: TSTEP_CONSTANT
      public :: TSTEP_INSTANT
      public :: TSTEP_AVG
      public :: TSTEP_ACCUM
      public :: TSTEP_MAX
      public :: TSTEP_MIN
      public :: TSTEP_DIFF
      public :: TSTEP_RMS
      public :: TSTEP_SD
      public :: TSTEP_COV
      public :: TSTEP_RATIO
      public :: TSTEP_RANGE
      public :: TSTEP_INSTANT2
      public :: TSTEP_INSTANT3
      public :: TAXIS_ABSOLUTE
      public :: TAXIS_RELATIVE
      public :: TAXIS_FORECAST
      public :: TUNIT_SECOND
      public :: TUNIT_MINUTE
      public :: TUNIT_QUARTER
      public :: TUNIT_30MINUTES
      public :: TUNIT_HOUR
      public :: TUNIT_3HOURS
      public :: TUNIT_6HOURS
      public :: TUNIT_12HOURS
      public :: TUNIT_DAY
      public :: TUNIT_MONTH
      public :: TUNIT_YEAR
      public :: CALENDAR_STANDARD
      public :: CALENDAR_PROLEPTIC
      public :: CALENDAR_360DAYS
      public :: CALENDAR_365DAYS
      public :: CALENDAR_366DAYS
      public :: CALENDAR_NONE
      public :: CDI_UUID_SIZE

contains
      function cdiStringError(cdiErrno)
        integer(kind=c_int), value :: cdiErrno
        interface
          function cdiStringError_c(cdiErrno) bind(c,name='cdiStringError')
            import :: c_ptr,c_int,c_char
            integer(kind=c_int), value :: cdiErrno
            type(c_ptr) :: cdiStringError_c
          end function cdiStringError_c
        end interface
        character(len=1, kind=c_char), pointer :: cdiStringError(:)
        type(c_ptr) :: cptr
        integer :: slen(1)

        cptr = cdiStringError_c(cdiErrno)
        cdiStringError => null()
        slen(1) = int(strlen(cptr))
        call c_f_pointer(cptr, cdiStringError, slen)
      end function cdiStringError
      function cdiLibraryVersion()
        interface
          function cdiLibraryVersion_c() bind(c,name='cdiLibraryVersion')
            import :: c_ptr,c_char
            type(c_ptr) :: cdiLibraryVersion_c
          end function cdiLibraryVersion_c
        end interface
        character(len=1, kind=c_char), pointer :: cdiLibraryVersion(:)
        type(c_ptr) :: cptr
        integer :: slen(1)

        cptr = cdiLibraryVersion_c()
        cdiLibraryVersion => null()
        slen(1) = int(strlen(cptr))
        call c_f_pointer(cptr, cdiLibraryVersion, slen)
      end function cdiLibraryVersion
      function streamFilename(streamID)
        integer(kind=c_int), value :: streamID
        interface
          function streamFilename_c(streamID) bind(c,name='streamFilename')
            import :: c_ptr,c_int,c_char
            integer(kind=c_int), value :: streamID
            type(c_ptr) :: streamFilename_c
          end function streamFilename_c
        end interface
        character(len=1, kind=c_char), pointer :: streamFilename(:)
        type(c_ptr) :: cptr
        integer :: slen(1)

        cptr = streamFilename_c(streamID)
        streamFilename => null()
        slen(1) = int(strlen(cptr))
        call c_f_pointer(cptr, streamFilename, slen)
      end function streamFilename
      function streamFilesuffix(filetype)
        integer(kind=c_int), value :: filetype
        interface
          function streamFilesuffix_c(filetype) bind(c,name='streamFilesuffix')
            import :: c_ptr,c_int,c_char
            integer(kind=c_int), value :: filetype
            type(c_ptr) :: streamFilesuffix_c
          end function streamFilesuffix_c
        end interface
        character(len=1, kind=c_char), pointer :: streamFilesuffix(:)
        type(c_ptr) :: cptr
        integer :: slen(1)

        cptr = streamFilesuffix_c(filetype)
        streamFilesuffix => null()
        slen(1) = int(strlen(cptr))
        call c_f_pointer(cptr, streamFilesuffix, slen)
      end function streamFilesuffix
      function gridNamePtr(gridtype)
        integer(kind=c_int), value :: gridtype
        interface
          function gridNamePtr_c(gridtype) bind(c,name='gridNamePtr')
            import :: c_ptr,c_int,c_char
            integer(kind=c_int), value :: gridtype
            type(c_ptr) :: gridNamePtr_c
          end function gridNamePtr_c
        end interface
        character(len=1, kind=c_char), pointer :: gridNamePtr(:)
        type(c_ptr) :: cptr
        integer :: slen(1)

        cptr = gridNamePtr_c(gridtype)
        gridNamePtr => null()
        slen(1) = int(strlen(cptr))
        call c_f_pointer(cptr, gridNamePtr, slen)
      end function gridNamePtr
      function tunitNamePtr(tunitID)
        integer(kind=c_int), value :: tunitID
        interface
          function tunitNamePtr_c(tunitID) bind(c,name='tunitNamePtr')
            import :: c_ptr,c_int,c_char
            integer(kind=c_int), value :: tunitID
            type(c_ptr) :: tunitNamePtr_c
          end function tunitNamePtr_c
        end interface
        character(len=1, kind=c_char), pointer :: tunitNamePtr(:)
        type(c_ptr) :: cptr
        integer :: slen(1)

        cptr = tunitNamePtr_c(tunitID)
        tunitNamePtr => null()
        slen(1) = int(strlen(cptr))
        call c_f_pointer(cptr, tunitNamePtr, slen)
      end function tunitNamePtr
      function institutInqNamePtr(instID)
        integer(kind=c_int), value :: instID
        interface
          function institutInqNamePtr_c(instID) bind(c,name='institutInqNamePtr')
            import :: c_ptr,c_int,c_char
            integer(kind=c_int), value :: instID
            type(c_ptr) :: institutInqNamePtr_c
          end function institutInqNamePtr_c
        end interface
        character(len=1, kind=c_char), pointer :: institutInqNamePtr(:)
        type(c_ptr) :: cptr
        integer :: slen(1)

        cptr = institutInqNamePtr_c(instID)
        institutInqNamePtr => null()
        slen(1) = int(strlen(cptr))
        call c_f_pointer(cptr, institutInqNamePtr, slen)
      end function institutInqNamePtr
      function institutInqLongnamePtr(instID)
        integer(kind=c_int), value :: instID
        interface
          function institutInqLongnamePtr_c(instID) bind(c,name='institutInqLongnamePtr')
            import :: c_ptr,c_int,c_char
            integer(kind=c_int), value :: instID
            type(c_ptr) :: institutInqLongnamePtr_c
          end function institutInqLongnamePtr_c
        end interface
        character(len=1, kind=c_char), pointer :: institutInqLongnamePtr(:)
        type(c_ptr) :: cptr
        integer :: slen(1)

        cptr = institutInqLongnamePtr_c(instID)
        institutInqLongnamePtr => null()
        slen(1) = int(strlen(cptr))
        call c_f_pointer(cptr, institutInqLongnamePtr, slen)
      end function institutInqLongnamePtr
      function modelInqNamePtr(modelID)
        integer(kind=c_int), value :: modelID
        interface
          function modelInqNamePtr_c(modelID) bind(c,name='modelInqNamePtr')
            import :: c_ptr,c_int,c_char
            integer(kind=c_int), value :: modelID
            type(c_ptr) :: modelInqNamePtr_c
          end function modelInqNamePtr_c
        end interface
        character(len=1, kind=c_char), pointer :: modelInqNamePtr(:)
        type(c_ptr) :: cptr
        integer :: slen(1)

        cptr = modelInqNamePtr_c(modelID)
        modelInqNamePtr => null()
        slen(1) = int(strlen(cptr))
        call c_f_pointer(cptr, modelInqNamePtr, slen)
      end function modelInqNamePtr
      function tableInqNamePtr(tableID)
        integer(kind=c_int), value :: tableID
        interface
          function tableInqNamePtr_c(tableID) bind(c,name='tableInqNamePtr')
            import :: c_ptr,c_int,c_char
            integer(kind=c_int), value :: tableID
            type(c_ptr) :: tableInqNamePtr_c
          end function tableInqNamePtr_c
        end interface
        character(len=1, kind=c_char), pointer :: tableInqNamePtr(:)
        type(c_ptr) :: cptr
        integer :: slen(1)

        cptr = tableInqNamePtr_c(tableID)
        tableInqNamePtr => null()
        slen(1) = int(strlen(cptr))
        call c_f_pointer(cptr, tableInqNamePtr, slen)
      end function tableInqNamePtr
      function tableInqParNamePtr(tableID,parID)
        integer(kind=c_int), value :: tableID
        integer(kind=c_int), value :: parID
        interface
          function tableInqParNamePtr_c(tableID,parID) bind(c,name='tableInqParNamePtr')
            import :: c_ptr,c_int,c_char
            integer(kind=c_int), value :: tableID
            integer(kind=c_int), value :: parID
            type(c_ptr) :: tableInqParNamePtr_c
          end function tableInqParNamePtr_c
        end interface
        character(len=1, kind=c_char), pointer :: tableInqParNamePtr(:)
        type(c_ptr) :: cptr
        integer :: slen(1)

        cptr = tableInqParNamePtr_c(tableID,&
          parID)
        tableInqParNamePtr => null()
        slen(1) = int(strlen(cptr))
        call c_f_pointer(cptr, tableInqParNamePtr, slen)
      end function tableInqParNamePtr
      function tableInqParLongnamePtr(tableID,parID)
        integer(kind=c_int), value :: tableID
        integer(kind=c_int), value :: parID
        interface
          function tableInqParLongnamePtr_c(tableID,parID) bind(c,name='tableInqParLongnamePtr')
            import :: c_ptr,c_int,c_char
            integer(kind=c_int), value :: tableID
            integer(kind=c_int), value :: parID
            type(c_ptr) :: tableInqParLongnamePtr_c
          end function tableInqParLongnamePtr_c
        end interface
        character(len=1, kind=c_char), pointer :: tableInqParLongnamePtr(:)
        type(c_ptr) :: cptr
        integer :: slen(1)

        cptr = tableInqParLongnamePtr_c(tableID,&
          parID)
        tableInqParLongnamePtr => null()
        slen(1) = int(strlen(cptr))
        call c_f_pointer(cptr, tableInqParLongnamePtr, slen)
      end function tableInqParLongnamePtr
      function tableInqParUnitsPtr(tableID,parID)
        integer(kind=c_int), value :: tableID
        integer(kind=c_int), value :: parID
        interface
          function tableInqParUnitsPtr_c(tableID,parID) bind(c,name='tableInqParUnitsPtr')
            import :: c_ptr,c_int,c_char
            integer(kind=c_int), value :: tableID
            integer(kind=c_int), value :: parID
            type(c_ptr) :: tableInqParUnitsPtr_c
          end function tableInqParUnitsPtr_c
        end interface
        character(len=1, kind=c_char), pointer :: tableInqParUnitsPtr(:)
        type(c_ptr) :: cptr
        integer :: slen(1)

        cptr = tableInqParUnitsPtr_c(tableID,&
          parID)
        tableInqParUnitsPtr => null()
        slen(1) = int(strlen(cptr))
        call c_f_pointer(cptr, tableInqParUnitsPtr, slen)
      end function tableInqParUnitsPtr

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

    end subroutine ctrim

    function c_len(s) result(i)
      character(kind=c_char), intent(in) :: s(:)
      integer :: i
      do i = 1, size(s)
        if (s(i) == c_null_char) then
          exit
        end if
      end do
      i = i - 1
    end function


end module mo_cdi

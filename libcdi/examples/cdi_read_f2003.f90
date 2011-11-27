PROGRAM CDIREADF2003
  use iso_c_binding
  use mo_cdi

  IMPLICIT NONE

  INTEGER :: gsize, nlevel, nvars, code
  INTEGER :: vdate, vtime, nmiss, status, ilev
  INTEGER :: streamID, varID, gridID, zaxisID
  INTEGER :: tsID, vlistID, taxisID
  DOUBLE PRECISION, ALLOCATABLE :: field(:,:)
  CHARACTER(kind=c_char,len=256) :: name, longname, units, msg

  ! Open the dataset
  streamID = streamOpenRead(C_CHAR_"example.nc"//C_NULL_CHAR)
  IF ( streamID < 0 ) THEN
    PRINT *,  'Could not Read the file.'
    msg = cdiStringError(streamID)
    WRITE(0,*) msg
    STOP 1
  END IF

  ! Get the variable list of the dataset
  vlistID = streamInqVlist(streamID)

  nvars = vlistNvars(vlistID)

  DO varID = 0, nvars-1
    code = vlistInqVarCode(vlistID, varID)
    CALL vlistInqVarName(vlistID, varID, name)
    CALL vlistInqVarLongname(vlistID, varID, longname)
    CALL vlistInqVarUnits(vlistID, varID, units)

    CALL ctrim(name)
    CALL ctrim(longname)
    CALL ctrim(units)

    WRITE(*,*) 'Parameter: ', varID+1, code,' ',trim(name),' ', &
                trim(longname),' ',trim(units), ' |'

  END DO

  ! Get the Time axis form the variable list
  taxisID = vlistInqTaxis(vlistID)

  ! Loop over the time steps
  DO tsID = 0, 999999
    ! Read the time step
    status = streamInqTimestep(streamID, tsID)
    IF ( status == 0 ) exit

    ! Get the verification date and time
    vdate = taxisInqVdate(taxisID)
    vtime = taxisInqVtime(taxisID)

    WRITE(*,*) 'Timestep: ', tsID+1, vdate, vtime

    ! Read the variables at the current timestep
    DO varID = 0, nvars-1
      gridID = vlistInqVarGrid(vlistID, varID)
      gsize = gridInqSize(gridID)
      zaxisID = vlistInqVarZaxis(vlistID, varID)
      nlevel = zaxisInqSize(zaxisID)
      ALLOCATE(field(gsize, nlevel))
      CALL streamReadVar(streamID, varID, field, nmiss)
      DO ilev = 1, nlevel
        WRITE(*,*) '   var=', varID+1, ' level=', ilev, ':', &
                  MINVAL(field(:,ilev)), MAXVAL(field(:,ilev))
      END DO
      DEALLOCATE(field)
    END DO
  END DO

  ! Close the input stream
  CALL streamClose(streamID)

END PROGRAM CDIREADF2003

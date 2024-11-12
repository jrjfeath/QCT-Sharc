# QCT-Sharc
 A modification of Sharc-3.0.2 for QCT calculations using NO+Ar

This project was developed to run QCT calculations using sharc and my other repository iScatter created in collaboration with Chris Robertson. 
This project requires the outputs folder produced by iScatter to be in the same directory as QCT-main.py to run

Plot-Geometries.py takes an output.dat produced in a trajectory and plots each step for NO+Ar

## Installation:
1. Download the zip file of sharc from https://github.com/sharc-md/sharc/releases/ 
2. Extract and upload to relevant server
3. Overwrite source folder with this source folder (Only do this for sharc-3.0.2)
4. If you must use a newer version of sharc you can make the edits listed below to hopefully maintain functionality.
5. After files have been overwritten run: make
6. make install
7. Sharc should now be ready to run QCT
8. You will need to edit QCT-Main.py and point to your installation directory
9. If you change the PES you will need to edit the relevant function calls in QCT-Main.py, no edits should need to be made to sharc

The provided PESV2.f is a subroutine that allows data to be obtained from the PES. It needs to be compiled using:
- python -m numpy.f2py -c PESV2.f -m PESV2 --compiler=intelem --fcompiler=intelem

The following calls in sharc (driver.f90, interface.F90, main.F90, output.f90) have been removed to greatly increase speed:
- call write_list_line(u_lis,traj,ctrl)
- call write_geom(u_geo, traj, ctrl)
- call write_restart_traj(u_rest,ctrl,traj)

### Changes to: qm.f90
The following subroutine has been completely rewritten:
```
subroutine call_runqm(traj)
  use definitions
  use output
  implicit none
  type(trajectory_type) :: traj
  character(len=100) :: data_from_python
  integer :: io_status
  read(*,*,iostat=io_status) data_from_python
  if (io_status /= 0) then
  return
  endif
endsubroutine
```

All calls related to writting QM.in have been removed:
```
open(u_qm_qmin,file='QM/QM.in',status='replace',action='write')
call write_infos(traj,ctrl)
call write_tasks_third(traj,ctrl)
close(u_qm_qmin)
```

And replaced with communication to python handled by replacing File I/O from above with:
```
write(6,*) 'QM Please'
do i=1,ctrl%natom
  write(6,'(A2,3(F12.7,1X),3X,3(F12.7,1X))') &
  &traj%element_a(i),(au2a*traj%geom_ad(i,j),j=1,3),(traj%veloc_ad(i,j),j=1,3)
enddo
write(6,*) 'Done!'
```
Replacing runQM.sh with stdin and stdout communication between python and fortran resulted in 100x speed increase!

### Changes to: output.f90
The following subroutine has been completely rewritten:
```
subroutine write_dat(u, traj, ctrl)
  use definitions
  use matrix
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  integer :: u, i, j
  integer :: nstates, natom, stride, io_status
  character(len=100) :: data_from_python

  ! check if writing
  stride=ctrl%output_steps_stride(1)
  if (traj%step>=ctrl%output_steps_limits(2)) then
    stride=ctrl%output_steps_stride(2)
  endif
  if (traj%step>=ctrl%output_steps_limits(3)) then
    stride=ctrl%output_steps_stride(3)
  endif
  if (modulo(traj%step,stride)==0) then

    nstates=ctrl%nstates
    natom=ctrl%natom

    ! Send geometry to Python
    write(6,*) 'GEOM'
    do i=1,natom
      write(6,*) (traj%geom_ad(i,j),j=1,3)
    enddo
    ! Wait for Python to say it got the data
    read(*,*,iostat=io_status) data_from_python
    if (io_status /= 0) then
      return
    endif
    ! Send velocity to Python
    write(6,*) 'VELOC'
    do i=1,natom
      write(6,*) (traj%veloc_ad(i,j),j=1,3)
    enddo
    ! Wait for Python to say it got the data
    read(*,*,iostat=io_status) data_from_python
    if (io_status /= 0) then
      return
    endif

  endif  ! <-- end of stride check

endsubroutine
```

### Changes to: qm_out.f90
- The data is now passed from python via stdout rather than files. 
- All get_thing routines need to have the goto_flag call updated to include the line_index
- - call goto_flag(2,'get_dipoles', line_index) <- line_index stores where in the buffer that property is found
- matread has been replaced with pyread, d3vecread has been replaced by pyd3vecread
```
subroutine goto_flag(flag1, routine, line_index)
  implicit none
  character(len=*) :: routine
  integer :: flag1, flag, i, io
  integer, intent(out) :: line_index
  character(len=200) :: string

  ! Loop over data_buffer
  do i = 1, n_lines
    string = trim(data_buffer(i))

    ! Check if the line starts with '!'
    if (string(1:1) == '!') then
      ! Try to read the flag number from the string
      read(string(2:20), *, iostat=io) flag
      if (io == 0 .and. flag == flag1) then
        ! Found the matching flag, set line_index
        line_index = i + 1  ! Move to the next line after the flag
        return
      end if
    end if
  end do

  ! If we get here, the flag wasn't found
  write(0,*) 'Quantity not found in data buffer'
  write(0,*) 'Routine =', trim(routine)
  stop 1

end subroutine

subroutine goto_flag_nostop(flag1, stat)
  implicit none

  integer :: flag1, flag, stat, i, io
  character(len=200) :: string

  ! Loop through the data_buffer
  do i = 1, n_lines
    string = trim(data_buffer(i))

    ! Check if the line starts with '!'
    if (string(1:1) == '!') then
      ! Try to read the flag number from the string
      read(string(2:20), *, iostat=io) flag
      if (io == 0 .and. flag == flag1) then
        ! Flag found, set stat to 0 and return
        stat = 0
        return
      end if
    end if
  end do

  ! Flag was not found, set stat to -1
  stat = -1
  return

end subroutine

subroutine pyzread(n, A, line_index, title)
  implicit none
  integer, intent(in) :: n
  integer, intent(inout) :: line_index
  complex*16, intent(out) :: A(n, n)
  character(len=8000), intent(out) :: title
  integer :: i, j, io
  real*8 :: line(2*n)
  character(len=200) :: string

  ! Read the title line from data_buffer
  title = trim(data_buffer(line_index))
  line_index = line_index + 1

  ! Read the matrix data from data_buffer
  do i = 1, n
    string = trim(data_buffer(line_index))
    line_index = line_index + 1
    read(string, *, iostat=io) (line(j), j = 1, 2*n)
    if (io /= 0) then
      write(*,*) 'Could not read matrix'
      write(*,*) 'routine=zread(), n=', n
      write(*,*) 'title=', trim(title)
      stop 1
    end if
    do j = 1, n
      A(i, j) = dcmplx(line(2*j-1), line(2*j))
    end do
  end do

end subroutine

subroutine pyd3vecread(n, c, line_index, title)
  implicit none
  ! parameters
  integer, intent(in) :: n
  real*8, intent(out) :: c(n, 3)
  integer, intent(inout) :: line_index  ! Current position in data_buffer
  character(len=8000), intent(out) :: title
  ! internal variables
  integer :: i, j, io
  character(len=200) :: string

  ! Read the title line from data_buffer
  title = trim(data_buffer(line_index))
  line_index = line_index + 1

  ! Read the 3-vectors from data_buffer
  do i = 1, n
    if (line_index > n_lines) then
      write(*,*) 'Error: Reached end of data_buffer while reading 3-vectors'
      stop 1
    end if
    string = trim(data_buffer(line_index))
    line_index = line_index + 1

    ! Read the 3-vector components from the line
    read(string, *, iostat=io) (c(i, j), j = 1, 3)
    if (io /= 0) then
      write(*,*) 'Could not read 3-vector'
      write(*,*) 'routine=d3vecread(), n=', n
      write(*,*) 'title=', trim(title)
      stop 1
    end if
  end do

end subroutine
```
# QCT-Sharc
 A modification of Sharc-3.0.2 for QCT calculations using NO+Ar

This project was developed to run QCT calculations using sharc and my other repository iScatter created in collaboration with Chris Robertson. 

Plot-Geometries.py takes an output.dat produced in a trajectory and plots each step for NO+Ar

The provided PESV2.f is a subroutine that allows data to be obtained from the PES. It needs to be compiled using:
- python -m numpy.f2py -c PESV2.f -m PESV2 --compiler=intelem --fcompiler=intelem

Installation:
1. Download the zip file of sharc from https://github.com/sharc-md/sharc/releases/ 
2. Extract and upload to relevant server
3. Overwrite source folder with this source folder (Only do this for sharc-3.0.2)
4. If you must use a newer version of sharc you can make the edits listed below to hopefully maintain functionality.
5. After files have been overwritten run: make
6. make install
7. Sharc should now be ready to run QCT
8. You will need to edit QCT-Main.py and point to your installation directory
9. If you change the PES you will need to edit the relevant function calls in QCT-Main.py, no edits should need to be made to sharc

The following calls in sharc (driver.f90, interface.F90, main.F90, output.f90) have been removed to greatly increase speed:
- call write_list_line(u_lis,traj,ctrl)
- call write_geom(u_geo, traj, ctrl)
- call write_restart_traj(u_rest,ctrl,traj)

The following changes have been made to qm.f90:
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

The following changes have been made to output.f90:
    The following subroutine has been completely rewritten:
    ```
    subroutine write_dat(u, traj, ctrl)
        use definitions
        use matrix
        implicit none
        type(trajectory_type) :: traj
        type(ctrl_type) :: ctrl
        integer :: u, i, j
        character(8000) :: string

        integer :: nstates, natom, stride

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

            if (ctrl%integrator==2) then 
            write(u,'(A)') '! 0 Step'
            write(u,'(I12)') traj%step
            else if (ctrl%integrator==1 .or. ctrl%integrator==0) then 
            write(u,'(A)') '! 0 Step'
            write(u,'(I12,F12.6,F12.6)') traj%step, traj%microtime*au2fs, ctrl%dtstep*au2fs
            endif

            call vec3write(natom, traj%geom_ad, u, '! 11 Geometry in a.u.','E21.13e3')
            call vec3write(natom, traj%veloc_ad, u, '! 12 Velocities in a.u.','E21.13e3')

        endif  ! <-- end of stride check
    endsubroutine
    ```
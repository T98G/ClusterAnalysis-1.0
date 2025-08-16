program clusterAnalysis
  
  implicit none
        
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Define some variables !!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer :: i, j, k
  integer :: status, atom_count, frame_count
  integer :: residue_count, temp
  integer :: start, finish, clock_rate, clock_max
  integer :: max_neighbors

  integer, allocatable :: adj_list(:, :)
  integer, allocatable :: neighbors_list(:, :, :)
  integer, allocatable :: cluster_counts(:)
  integer, allocatable :: main_cluster_sizes(:)

  real :: hour, minute, second
  real, allocatable :: cluster_sizes(:)
  real, allocatable :: neighbors_probability(:)
 
  double precision :: x, y, z
  double precision :: radius, distance
  double precision :: box_dimensions(3)
  double precision, allocatable :: coordinates(:, :)

  character(len=30), dimension(30) :: args
  character(len=3) :: selection
  character(len=100) :: neighbors_file
  character(len=100) :: pdb, counts_file, hist_file
  character(len=100) :: line, sizes_file, npd_file
  
  logical, allocatable :: visited(:)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initialize some variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pdb = ""
  selection = ""
  counts_file = ""
  hist_file = ""
  neighbors_file = ""
  sizes_file = ""
  npd_file = ""
  radius = 0
  max_neighbors = 0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!! Get and treat command-line arguments !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i = 1, 30        
     call getarg(i, args(i))  
  end do
  
  do i = 1, 30
    
    if (args(i) == "-help" .or. args(i) == "-h" .or. args(i) == "--h") then      
      call usage()
      stop
    end if
    
    if (args(i) == "-pdb") then
       pdb = trim(args(i + 1))
    end if
    
    if (args(i) == "-radius") then
      read(args(i + 1), *) radius
    end if

    if (args(i) == "-selection") then
      selection = trim(args(i + 1)) 
    end if

    if (args(i) == "-sizes") then
      sizes_file = trim(args(i + 1))
    end if

    if (args(i) == "-counts") then
      counts_file = trim(args(i + 1))      
    end if

    if (args(i) == "-hist") then
      hist_file = trim(args(i + 1))
    end if

    if (args(i) == "-n") then
      read(args(i + 1), *) max_neighbors
    end if

    if (args(i) == "-neighbors") then
       neighbors_file = trim(args(i + 1))
    end if

    if (args(i) == "-npd") then
      npd_file = trim(args(i + 1))
    end if
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Check Comand Line arguments !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if (pdb == "") then    
    write(*, *) "PDB file was not provided! Stoping."
    stop
  end if
  
  if (radius == 0) then    
    write(*, *) "The cutoff radius cannot be 0! Stoping."
    stop
  end if

  if (hist_file == "") hist_file = "cluster_hist.xvg"

  if (counts_file == "") counts_file = "cluster_counts.xvg"
  
  if (max_neighbors == 0) then
    
    write(*, *) "Setting Maximum neighbors to 15"
    
    max_neighbors = 15

  end if

  if (neighbors_file == "") neighbors_file = "neighbors.dat"
  
  if (sizes_file == "") sizes_file = "size.xvg"

  if (npd_file == "") npd_file = "neighbors_prob.xvg"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initialize counts !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  atom_count = 0
  frame_count = 0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Open file for reading !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  open(unit=1, file=pdb, status='old', action='read')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Start Clock !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call system_clock(start, clock_rate, clock_max)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!! Determine File Size !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do     
    read(1, '(A)', iostat=status) line     
    if (status /= 0) exit  ! Exit loop if end of file or an error

    if (frame_count == 0 .and. line(1:4) == "ATOM") then        

      if (index(line(12:16), selection) /= 0) then
        atom_count = atom_count + 1
        read(line(23:26), '(I4)') residue_count
      end if
    
    end if
    
    if (line(1:3) == 'END') then
      frame_count = frame_count + 1
    end if

  end do

  !Make sure there is at least one frame   
  frame_count = max(frame_count, 1)
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Allocate memory !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(coordinates(3, atom_count))
  allocate(neighbors_list(atom_count, atom_count, frame_count))
  allocate(adj_list(max_neighbors, atom_count))
  allocate(visited(atom_count))
  allocate(cluster_counts(frame_count))
  allocate(cluster_sizes(atom_count))
  allocate(main_cluster_sizes(frame_count))
  allocate(neighbors_probability(max_neighbors))
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initialize Variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  coordinates = 0
  box_dimensions = 0
  neighbors_list = 0
  adj_list = 0
  cluster_counts = 0
  cluster_sizes = 0
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Open file for Writting !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(unit=9, file=neighbors_file, status="new")

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Main loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  rewind(1)
  
  write(*, *)
  write(*, *)
  write(*, *) "Calculating"

  do i = 1, frame_count
      
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!! Loop Through the file and read lines !!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    coordinates = 0
    box_dimensions = 0
    adj_list = 0

    j = 1       
    !! Defining J outside the loop initialization because there are header lines
    !! Which are not atoms so I only want to increment j if the line is an atom

    do while (.true.)       
      read(1, '(A)', iostat=status) line
      if (status /= 0) exit
      if (line(1:3) == "END") exit

      if (line(1:4) == "ATOM" .and. index(line(12:16), selection) /= 0) then            
        
        read(line(31:54), *) x, y, z
        coordinates(1, j) = x
        coordinates(2, j) = y
        coordinates(3, j) = z
        j = j + 1 
      end if

      if (line(1:6) == 'CRYST1') then    
        read(line(7:15), *) box_dimensions(1)
        read(line(16:24), *) box_dimensions(2)
        read(line(25:33), *) box_dimensions(3)    
      end if
        
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!! Create Adjacency list !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do j = 1, atom_count - 1

      do k = j + 1, atom_count

        distance = Calculate_distance(coordinates(:, j), coordinates(:, k), box_dimensions)     
        
        if (distance < radius) then
          ! Add k to the adjacency list of j
          call AddToAdjList(adj_list(:, j), k, max_neighbors)
          ! Add j to the adjacency list of k
          call AddToAdjList(adj_list(:, k), j, max_neighbors)                     
        end if

      end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!! Write adjacency list to file !!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    write(9,*) 'Array of neighbors, frame:', i
    do j = 1, atom_count
      write(9,*) (adj_list(k,j),k=1,max_neighbors)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! Make Number of neighbors Histogram !!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    do j = 1, atom_count

      temp = (count(adj_list(:, j) /= 0))

      if (temp /=0) then 
        neighbors_probability(temp) = neighbors_probability(temp) + 1
      end if

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Find Neighbors !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Set visited array to false 
    visited = .false.
    
    do j = 1, atom_count      
      call FindIndirectNeighbors(j, adj_list, neighbors_list(:, j, i), max_neighbors, visited)    
    end do
       
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! Count Unique Clusters !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !We start searching recursively the neighbors of atom j 
    !Considering that I'm only populating the neighbors list for atom j, all the neighbors lists for
    !any atom in the same cluster of atom j, will contain only zeros
    !so if a list contain a non zero value it's for sure another cluster
        
    do j = 1, atom_count
      
      if (count(neighbors_list(:, j, i) /= 0) > 0) cluster_counts(i) = cluster_counts(i) + 1

    end do 
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!! Print Progress Bar !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call progressBar(real(i) / real(frame_count))

  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Close some files !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  close(1)
  close(9)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!! Calculate Main Cluster Size !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  open(unit=12, file=sizes_file, status="new")
  
  write(12, *) '@    title "Main Cluster Size"'
  write(12, *) '@    xaxis  label "Size"'
  write(12, *) '@    yaxis  label "Frame Number"'
  write(12, *) '@    legend off'
    
  do i = 1, frame_count
    do j = 1, atom_count

      cluster_sizes(j) = count(neighbors_list(:, j, i) /= 0)

    end do

      write(12, '(I9, F10.5)') i, maxval(cluster_sizes)

  end do
  
  close(12)

  cluster_sizes(:) = 0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!! Write neighbors count histogram to file !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  neighbors_probability(:) = neighbors_probability(:) / sum(neighbors_probability(:))

  open(unit=13, file=npd_file, status="new")
  
  write(13, *) '@    title "Neighbors Count Probability"'
  write(13, *) '@    xaxis  label "Number of Neighbors"'
  write(13, *) '@    yaxis  label "Probability"'
  write(13, *) '@    legend off'

  do i = 1, max_neighbors
    
    write(13, "(I9, F10.5)") i, neighbors_probability(i)

  end do
  
  close(13)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Generate Histogram !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  temp = 0

  do i = 1, frame_count
    do j = 1, atom_count
      
      temp = count(neighbors_list(:, j, i) /= 0)  
      
      if (temp > 0) then
        cluster_sizes(temp) = cluster_sizes(temp) + 1
      end if

    end do
  end do

  cluster_sizes = cluster_sizes / sum(cluster_sizes)
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Save Results to Files !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  ! Save results to two xvg files for XMGrace
  
  ! Save histogram file

  open(unit=10, file=hist_file, status="new")

  write(10, *) '@    title "Cluster Size Distribution"'
  write(10, *) '@    xaxis  label "Cluster Size"'
  write(10, *) '@    yaxis  label "Probability"'
  write(10, *) '@    legend off'
    
  do i = 1, atom_count
      
    write(10, '(I9, 1X, F10.5)') i, cluster_sizes(i)

  end do

  close(10)
 
  !Save Cluster Counts

  open(unit=15, file=counts_file, status="new")
    
  write(15, *) '@    title "Cluster Counts"'
  write(15, *) '@    xaxis  label "Frame Number"'    
  write(15, *) '@    yaxis  label "Number of Clusters"'
  write(15, *) '@    legend off'
    
  do i = 1, frame_count
    write(15, '(I9, 1X, I9)') i, cluster_counts(i) 
  end do

  close(15)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Deallocate Memory !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  deallocate(coordinates)
  deallocate(visited)
  deallocate(adj_list)
  deallocate(neighbors_list)
  deallocate(cluster_counts)
  deallocate(cluster_sizes) 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Stop clock !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  call system_clock(finish, clock_rate, clock_max)

  hour = real(finish - start) / (3600 * real(clock_rate))
  minute = (hour - int(hour)) * 60
  second = (minute - int(minute)) * 60

  write(*, *)
  write(*, *)

  write(*, 1) " Processing time : ", int(hour), " h", int(minute), " min", int(second), " sec"
              1 format(a19, i6.3, a2, i6.3, a4, i6.3, a4)

  write(*, *)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Define Subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains

  subroutine usage()
    
    write(*, *) ""
    write(*, *) "Welcome to cluster analysis" 
    write(*, *) "" 
    write(*, *) "This program calculates a few things:"
    write(*, *) ""
    write(*, *) " * The number of clusters per frame."
    write(*, *) " * A histogram of the cluster sizes probability."
    write(*, *) " * A file with a list of neighbors per atom"
    write(*, *) " * The cluster size of the largest cluster as a function of time"
    write(*, *) ""
    write(*, *) ""
    write(*, *) ""
    write(*, *) "Usage:"
    write(*, *) ""
    write(*, *) "-pdb       The input pdb trajectory file"
    write(*, *) "-radius    The radius in Angstroms for clustering"
    write(*, *) "-selection The atom name selection for clustering"
    write(*, *) ""
    write(*, *) "Optional flags:"
    write(*, *) ""
    write(*, *) "-counts    The name for the cluster counts output file"
    write(*, *) "-sizes     The name for the main cluster size output file"
    write(*, *) "-hist      The name for the histogram output file"
    write(*, *) "-n         The Maximum number of neighbors" 
    write(*, *) "-neighbors A .dat output file with the list of neighbors for each atom"
    write(*, *) ""
    write(*, *) ""

  end subroutine usage


  subroutine progressBar(progress)
    
    implicit none
    real, intent(in) :: progress
    integer :: PBWIDTH
    character(len=60) :: PBSTR
    character(len=60) :: CTPBSTR
    integer :: lpad
    integer :: rpad

    PBWIDTH = 60
    PBSTR = "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"

    ! Ensure progress is between 0.0 and 1.0
    if (progress < 0.0) then
        lpad = 0
    else if (progress > 1.0) then
        lpad = PBWIDTH
    else
        lpad = int(progress * PBWIDTH)
    end if

    rpad = PBWIDTH - lpad

    ! Ensure lpad + rpad is equal to PBWIDTH
    if (lpad + rpad /= PBWIDTH) then
        lpad = PBWIDTH - rpad
    end if

    ! Construct the current progress bar string
    CTPBSTR = PBSTR(1:lpad)
    if (rpad > 0) then
        CTPBSTR(lpad + 1:lpad + rpad) = ' '
    end if

    ! Output the progress bar
    write(*,'(a, a, f5.1, a2, a1, a1)', advance="no") &
            char(13), ' progress: ', progress*100, '%'
    write(*, "(A)", advance="no") "["
    write(*, "(A)", advance="no")  CTPBSTR
    write(*, "(A)", advance="no") "]"
    flush(6)

  end subroutine progressBar

  function calculate_distance(coordinates_i, coordinates_j, box_dimensions) result(distance)
   
    implicit none
    
    double precision :: coordinates_i(3), coordinates_j(3)
    double precision :: box_dimensions(3) 
    double precision :: dist_pbc(3)
    double precision :: distance
    integer :: i
    
    ! Compute the distance
    
    do i = 1, 3
      dist_pbc(i) = abs(coordinates_i(i) - coordinates_j(i))
      if (dist_pbc(i) > box_dimensions(i) / 2) then
        dist_pbc(i) = abs(box_dimensions(i) - dist_pbc(i))
      end if
    end do
    
    distance = sqrt(dist_pbc(1) ** 2 + dist_pbc(2) ** 2 + dist_pbc(3) ** 2)
     
  end function calculate_distance


  subroutine AddToAdjList(adj_list_row, node, max_neighbors)
    
    implicit none

    integer, intent(inout) :: adj_list_row(:)
    integer, intent(in) :: node
    integer, intent(in) :: max_neighbors
    integer :: i

    do i = 1, max_neighbors
      
      if (adj_list_row(i) == node) return
      
      if (adj_list_row(i) == 0) then
        
         adj_list_row(i) = node
         return
      
      end if
    end do
  
  write(*, *) "Adjacency list is full, cannot add more atoms"  

  end subroutine AddToAdjList
  

  subroutine AddToNeighborsList(list_row, neighbors_list_frame_row, max_neighbors)

    implicit none
    
    integer, intent(in) :: list_row(:)
    integer, intent(in) :: max_neighbors
    integer, intent(inout) :: neighbors_list_frame_row(:)
    integer :: i, j

    do i = 1, max_neighbors  
        
        if (.not. any(neighbors_list_frame_row == list_row(i))) then

          do j = 1, size(neighbors_list_frame_row)
            
            if (neighbors_list_frame_row(j) == 0) then

              neighbors_list_frame_row(j) = list_row(i)

              exit
            
            end if        

          end do
        end if
    end do

  end subroutine AddToNeighborsList


  recursive subroutine FindIndirectNeighbors(start, adj_list_frame, neighbors_list_frame_row, max_neighbors, visited)
    
    !What this subroutine does
    !1. Use an atom as start
    !2. If start is zero or has been visited exit
    !3. Add all values in the adjecency list of the start atom to the neighbor list
    !4. Loop through the adjecency list for the start atom
    !5. If the value in the adjecy list is not zero, use it as start for recursive neighbor searching

    implicit none

    integer, intent(in) :: start
    integer, intent(in) :: adj_list_frame(:, :)
    integer, intent(in) :: max_neighbors
    integer, intent(inout) :: neighbors_list_frame_row(:)
    logical, intent(inout) :: visited(size(adj_list_frame, 2))
    integer :: i
    
    ! Base case: return if start is zero or already visited
    if (start == 0 .or. visited(start)) return
    
    ! Mark this node as visited
    visited(start) = .true.

    ! Add the current adjacency list to the neighbors list
    call AddToNeighborsList(adj_list_frame(:, start), neighbors_list_frame_row, max_neighbors)

    ! Loop through the adjacency list for the current start atom
    do i = 1, max_neighbors
      if (adj_list_frame(i, start) /= 0) then
        
        ! Recursively find indirect neighbors
        call FindIndirectNeighbors(adj_list_frame(i, start), adj_list_frame, neighbors_list_frame_row, max_neighbors, visited)
      
      end if
    end do

  end subroutine FindIndirectNeighbors
  
end program clusterAnalysis

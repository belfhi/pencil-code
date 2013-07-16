; +
; NAME:
;       PC_READ_SUBVOL_RAW
;
; PURPOSE:
;       Read sub-volumes from var.dat, or other VAR files in an efficient way!
;
;       Returns one array from a snapshot (var) file generated by a
;       Pencil Code run, and another array with the variable labels.
;       If you need to be efficient, please use 'pc_collect.x' to combine
;       distributed varfiles before reading them in IDL.
;       This routine can also read reduced datasets from 'pc_reduce.x'.
;
; CATEGORY:
;       Pencil Code, File I/O
;
; CALLING SEQUENCE:
;       pc_read_subvol_raw, object=object, varfile=varfile, tags=tags,       $
;           datadir=datadir, var_list=var_list, varcontent=varcontent,       $
;           start_param=start_param, run_param=run_param, allprocs=allprocs, $
;           xs=xs, xe=xe, ys=ys, ye=ye, zs=zs, ze=ze, /addghosts, /trimall,  $
;           dim=dim, sub_dim=sub_dim, grid=grid, sub_grid=sub_grid,          $
;           time=time, name=name, /quiet, /swap_endian, /f77, /reduced
;
; KEYWORD PARAMETERS:
;    datadir: Specifies the root data directory. Default: './data'.  [string]
;    varfile: Name of the snapshot. Default: 'var.dat'.              [string]
;       time: Timestamp of the snapshot.                             [float]
;       name: Name to be used for the generated sub_grid structure.  [string]
;   allprocs: Load distributed (0) or collective (1 or 2) varfiles.  [integer]
;   /reduced: Load previously reduced collective varfiles (implies allprocs=1).
;        dim: Dimension structure of the global 3D-setup.            [struct]
;    sub_dim: Returns a dimension structure of the loaded sub-volume.[struct]
;       grid: Grid structure of the global 3D-setup.                 [struct]
;   sub_grid: Returns a grid structure of the loaded sub-volume.     [struct]
;
;      xs/xe: starting/ending x-coordinate of the sub-volume.        [integer]
;      ys/ye: starting/ending y-coordinate of the sub-volume.        [integer]
;      zs/ze: starting/ending z-coordinate of the sub-volume.        [integer]
;             Default values are the first/last physical grid points [l1:l2,m1:m2,n1:n2].
;
;     object: Array with the loaded data is retured.                 [4D-array]
;       tags: Structure of tag names and their indices withhin data. [structure]
;   var_list: Array of varcontent idlvars to read (default = all).   [string(*)]
;
; /addghosts: Adds ghost layers to the given x/y/z starting/ending coordinates.
;   /trimall: Remove ghost layers from the returned data, dim and grid.
;     /quiet: Suppress any information messages and summary statistics.
;
; EXAMPLES:
;
; * Load a sub-volume and display it in a GUI:
;       pc_read_subvol_raw, obj=vars, tags=tags, var_list=["lnrho","uu"], xs=10, xe=14, ys=11, ye=49
;       cmp_cslice, { uz:vars[*,*,*,tags.uz], lnrho:vars[*,*,*,tags.lnrho] }
;
; * Load a sub-volume including ghost cells and compute some physical quantity:
;       pc_read_subvol_raw, obj=vars, tags=tags, sub_dim=dim, sub_grid=grid, xs=10, xe=14, ys=11, ye=49, /addghosts
;       HR_ohm = pc_get_quantity ('HR_ohm', vars, tags, dim=dim)
;       tvscl, HR_ohm[*,*,0]          ; Display lowest physical z-layer.
;       tvscl, HR_ohm[dim.nx-1,*,*]   ; Display rightmost physical yz-cut.
;       cslice, HR_ohm                ; Explore 3D sub-volume in a GUI.
;
; MODIFICATION HISTORY:
;       $Id$
;       Adapted from: pc_read_slice_raw.pro, 4th May 2012
;
;-
pro pc_read_subvol_raw, object=object, varfile=varfile, tags=tags, datadir=datadir, var_list=var_list, varcontent=varcontent, start_param=start_param, run_param=run_param, trimall=trimall, allprocs=allprocs, reduced=reduced, xs=xs, xe=xe, ys=ys, ye=ye, zs=zs, ze=ze, addghosts=addghosts, dim=dim, sub_dim=sub_dim, grid=grid, sub_grid=sub_grid, time=time, name=name, quiet=quiet, swap_endian=swap_endian, f77=f77

	; Use common block belonging to derivative routines etc. so they can be set up properly.
	common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
	common cdat_limits, l1, l2, m1, m2, n1, n2, nx, ny, nz
	common cdat_grid, dx_1, dy_1, dz_1, dx_tilde, dy_tilde, dz_tilde, lequidist, lperi, ldegenerated
	common pc_precision, zero, one
	common cdat_coords, coord_system

	; Default settings.
	default, allprocs, 0
	default, reduced, 0
	default, addghosts, 0
	default, swap_endian, 0
	if (keyword_set (name)) then name += "_" else name = "pc_read_subvol_raw_"
	if (keyword_set (reduced)) then allprocs = 1

	; Default data directory.
	if (not keyword_set (datadir)) then datadir = pc_get_datadir()

	; Name and path of varfile to read.
	if (not keyword_set (varfile)) then varfile = 'var.dat'

	; Automatically switch to allprocs, if necessary.
	if (not keyword_set (allprocs) and not file_test (datadir+'/proc0/'+varfile, /regular)) then allprocs = 1

	; Set f77 keyword according to allprocs.
	if (keyword_set (allprocs)) then default, f77, 0
	default, f77, 1

	; Get necessary dimensions quietly.
	if (n_elements (dim) eq 0) then pc_read_dim, object=dim, datadir=datadir, reduced=reduced, /quiet

	; ... and check pc_precision is set for all Pencil Code tools.
	pc_set_precision, dim=dim, quiet=quiet

	; Local shorthand for some parameters.
	precision = dim.precision
	nxgrid = dim.nxgrid
	nygrid = dim.nygrid
	nzgrid = dim.nzgrid
	nwgrid = nxgrid * nygrid * nzgrid
	mxgrid = dim.mxgrid
	mygrid = dim.mygrid
	mzgrid = dim.mzgrid
	mwgrid = mxgrid * mygrid * mzgrid
	nghostx = dim.nghostx
	nghosty = dim.nghosty
	nghostz = dim.nghostz

	if (addghosts) then begin
		default, xs, 0
		default, ys, 0
		default, zs, 0
		default, xe, xs + nxgrid - 1
		default, ye, ys + nygrid - 1
		default, ze, zs + nzgrid - 1
		xe = xe + 2 * nghostx
		ye = ye + 2 * nghosty
		ze = ze + 2 * nghostz
	end else begin
		default, xs, 0
		default, ys, 0
		default, zs, 0
		default, xe, mxgrid - 1
		default, ye, mygrid - 1
		default, ze, mzgrid - 1
	end
	xns = xs + nghostx
	xne = xe - nghostx
	yns = ys + nghosty
	yne = ye - nghosty
	zns = zs + nghostz
	zne = ze - nghostz
	gx_delta = xe - xs + 1
	gy_delta = ye - ys + 1
	gz_delta = ze - zs + 1

	if (any ([xs, xne-xns, ys, yne-yns, zs, zne-zns] lt 0) or any ([xe, ye, ze] ge [mxgrid, mygrid, mzgrid])) then $
		message, 'pc_read_subvol_raw: sub-volume indices are invalid.'

	; Get necessary parameters quietly.
	if (n_elements (start_param) eq 0) then $
		pc_read_param, object=start_param, dim=dim, datadir=datadir, /quiet
	if (n_elements (run_param) eq 0) then begin
		if (file_test (datadir+'/param2.nml')) then begin
			pc_read_param, object=run_param, /param2, dim=dim, datadir=datadir, /quiet
		endif else begin
			print, 'Could not find '+datadir+'/param2.nml'
		endelse
	endif

	if (n_elements (grid) eq 0) then $
		pc_read_grid, object=grid, dim=dim, param=start_param, datadir=datadir, allprocs=allprocs, reduced=reduced, /quiet

	; Set the coordinate system.
	coord_system = start_param.coord_system

	; Read local dimensions.
	ipx_start = 0
	ipy_start = 0
	ipz_start = 0
	if (allprocs eq 1) then begin
		procdim = dim
		ipx_end = 0
		ipy_end = 0
		ipz_end = 0
	end else begin
		pc_read_dim, object=procdim, proc=0, datadir=datadir, /quiet
		if (allprocs eq 2) then begin
			ipx_end = 0
			ipy_end = 0
			procdim.nx = procdim.nxgrid
			procdim.ny = procdim.nygrid
			procdim.mx = procdim.mxgrid
			procdim.my = procdim.mygrid
			procdim.mw = procdim.mx * procdim.my * procdim.mz
		end else begin
			ipx_start = xs / procdim.nx
			ipy_start = ys / procdim.ny
			ipx_end   = (xe - 2 * nghostx) / procdim.nx
			ipy_end   = (ye - 2 * nghosty) / procdim.ny
		end
		ipz_start = zs / procdim.nz
		ipz_end   = (ze - 2 * nghostz) / procdim.nz
	end

	; Generate dim structure of the sub-volume.
	sub_dim = dim
	sub_dim.mx = gx_delta
	sub_dim.my = gy_delta
	sub_dim.mz = gz_delta
	sub_dim.nx = gx_delta - 2*nghostx
	sub_dim.ny = gy_delta - 2*nghosty
	sub_dim.nz = gz_delta - 2*nghostz
	sub_dim.mxgrid = sub_dim.mx
	sub_dim.mygrid = sub_dim.my
	sub_dim.mzgrid = sub_dim.mz
	sub_dim.nxgrid = sub_dim.nx
	sub_dim.nygrid = sub_dim.ny
	sub_dim.nzgrid = sub_dim.nz
	sub_dim.mw = sub_dim.mx * sub_dim.my * sub_dim.mz
	sub_dim.l1 = nghostx
	sub_dim.m1 = nghosty
	sub_dim.n1 = nghostz
	sub_dim.l2 = sub_dim.mx - nghostx - 1
	sub_dim.m2 = sub_dim.my - nghosty - 1
	sub_dim.n2 = sub_dim.mz - nghostz - 1

	; Local shorthand for some parameters.
	nx = sub_dim.nx
	ny = sub_dim.ny
	nz = sub_dim.nz
	nw = nx * ny * nz
	mx = sub_dim.mx
	my = sub_dim.my
	mz = sub_dim.mz
	mw = mx * my * mz
	l1 = sub_dim.l1
	l2 = sub_dim.l2
	m1 = sub_dim.m1
	m2 = sub_dim.m2
	n1 = sub_dim.n1
	n2 = sub_dim.n2

	; Read meta data and set up variable/tag lists.
	if (n_elements (varcontent) eq 0) then $
		varcontent = pc_varcontent(datadir=datadir,dim=dim,param=start_param,quiet=quiet)
	totalvars = (size(varcontent))[1]
	if (n_elements (var_list) eq 0) then begin
		var_list = varcontent[*].idlvar
		var_list = var_list[where (var_list ne "dummy")]
	endif

	; Display information about the files contents.
	content = ''
	for iv=0L, totalvars-1L do begin
		content += ', '+varcontent[iv].variable
		; For vector quantities skip the required number of elements of the f array.
		iv += varcontent[iv].skip
	endfor
	content = strmid (content, 2)

	tags = { time:0.0d0 }
	read_content = ''
	indices = [ -1 ]
	num_read = 0
	num = n_elements (var_list)
	for ov=0L, num-1L do begin
		tag = var_list[ov]
		iv = where (varcontent[*].idlvar eq tag)
		if (iv ge 0) then begin
			if (tag eq "uu") then begin
				tags = create_struct (tags, "uu", [num_read, num_read+1, num_read+2])
				tags = create_struct (tags, "ux", num_read, "uy", num_read+1, "uz", num_read+2)
				indices = [ indices, iv, iv+1, iv+2 ]
				num_read += 3
			endif else if (tag eq "aa") then begin
				tags = create_struct (tags, "aa", [num_read, num_read+1, num_read+2])
				tags = create_struct (tags, "ax", num_read, "ay", num_read+1, "az", num_read+2)
				indices = [ indices, iv, iv+1, iv+2 ]
				num_read += 3
			endif else begin
				tags = create_struct (tags, tag, num_read)
				indices = [ indices, iv ]
				num_read++
			endelse
			read_content += ', '+varcontent[iv].variable
		endif
	endfor
	read_content = strmid (read_content, 2)
	if (not keyword_set (quiet)) then begin
		print, ''
		print, 'The file '+varfile+' contains: ', content
		if (strlen (read_content) lt strlen (content)) then print, 'Will read only: ', read_content
		print, ''
		print, 'The sub-volume dimension is ', sub_dim.mx, sub_dim.my, sub_dim.mz
		print, ''
	endif
	if (f77 eq 0) then markers = 0 else markers = 1

	if (num_read le 0) then begin
		if (not keyword_set (quiet)) then message, 'WARNING: nothing to read!'
	end else begin
		indices = indices[where (indices ge 0)]

		; Initialize output buffer.
		if (precision eq 'D') then begin
			bytes = 8
			object = dblarr (gx_delta, gy_delta, gz_delta, num_read)
		end else begin
			bytes = 4
			object = fltarr (gx_delta, gy_delta, gz_delta, num_read)
		end
	end

	; Iterate over processors.
	for ipz = ipz_start, ipz_end do begin
		if (num_read le 0) then continue
		for ipy = ipy_start, ipy_end do begin
			for ipx = ipx_start, ipx_end do begin
				; Set sub-volume parameters.
				if (ipx eq ipx_start) then px_start = xs - ipx * procdim.nx else px_start = nghostx
				if (ipy eq ipy_start) then py_start = ys - ipy * procdim.ny else py_start = nghosty
				if (ipz eq ipz_start) then pz_start = zs - ipz * procdim.nz else pz_start = nghostz
				if (ipx eq ipx_end) then px_end = xe - ipx * procdim.nx else px_end = procdim.mx - 1 - nghostx
				if (ipy eq ipy_end) then py_end = ye - ipy * procdim.ny else py_end = procdim.my - 1 - nghosty
				if (ipz eq ipz_end) then pz_end = ze - ipz * procdim.nz else pz_end = procdim.mz - 1 - nghostz
				px_delta = px_end - px_start + 1
				py_delta = py_end - py_start + 1
				pz_delta = pz_end - pz_start + 1

				; Initialize read buffer.
				if (precision eq 'D') then buffer = dblarr (px_delta) else buffer = fltarr (px_delta)

				; Initialize processor specific parameters.
				iproc = ipx + ipy*dim.nprocx + ipz*dim.nprocx*dim.nprocy
				if (ipx eq ipx_start) then x_off = 0 else x_off = nghostx + ipx * procdim.nx - xs
				if (ipy eq ipy_start) then y_off = 0 else y_off = nghosty + ipy * procdim.ny - ys
				if (ipz eq ipz_start) then z_off = 0 else z_off = nghostz + ipz * procdim.nz - zs
				if (allprocs eq 1) then begin
					if (keyword_set (reduced)) then procdir = 'reduced' else procdir = 'allprocs'
				end else begin
					procdir = 'proc' + strtrim (iproc, 2)
				end

				; Build the full path and filename.
				filename = datadir+'/'+procdir+'/'+varfile

				; Check for existence and read the data.
				if (not file_test(filename)) then message, 'ERROR: File not found "'+filename+'"'

				; Open a varfile and read some data!
				openr, lun, filename, swap_endian=swap_endian, /get_lun
				mx = long64 (procdim.mx)
				mxy = mx * procdim.my
				mxyz = mxy * procdim.mz
				for pos = 0, num_read-1 do begin
					pa = indices[pos]
					for pz = pz_start, pz_end do begin
						for py = py_start, py_end do begin
							point_lun, lun, bytes * (px_start + py*mx + pz*mxy + pa*mxyz) + long64 (markers*4)
							readu, lun, buffer
							object[x_off:x_off+px_delta-1,y_off+py-py_start,z_off+pz-pz_start,pos] = buffer
						endfor
					endfor
				endfor
				close, lun
				free_lun, lun
			end
		end
	end

	; Tidy memory a little.
	undefine, buffer

	; Read timestamp.
	pc_read_var_time, time=t, varfile=varfile, datadir=datadir, allprocs=allprocs, reduced=reduced, procdim=procdim, param=start_param, /quiet

	; Crop grid.
	x = grid.x[xs:xe]
	y = grid.y[ys:ye]
	z = grid.z[zs:ze]

	; Prepare for derivatives.
	dx = grid.dx
	dy = grid.dy
	dz = grid.dz
	dx_1 = grid.dx_1[xs:xe]
	dy_1 = grid.dy_1[ys:ye]
	dz_1 = grid.dz_1[zs:ze]
	dx_tilde = grid.dx_tilde[xs:xe]
	dy_tilde = grid.dy_tilde[ys:ye]
	dz_tilde = grid.dz_tilde[zs:ze]
	if (xe-xs eq mxgrid-1) then Lx = grid.Lx else Lx = total (1.0/grid.dx_1[xns:xne])
	if (ye-ys eq mygrid-1) then Ly = grid.Ly else Ly = total (1.0/grid.dy_1[yns:yne])
	if (ze-zs eq mzgrid-1) then Lz = grid.Lz else Lz = total (1.0/grid.dz_1[zns:zne])
	ldegenerated = [ xe-xs, ye-ys, ze-zs ] lt 2 * [ nghostx, nghosty, nghostz ]

	if (keyword_set (reduced)) then name += "reduced_"

	; Remove ghost zones if requested.
	if (keyword_set (trimall)) then begin
		if (num_read ge 1) then object = object[l1:l2,m1:m2,n1:n2,*]
		x = x[l1:l2]
		y = y[m1:m2]
		z = z[n1:n2]
		dx_1 = dx_1[l1:l2]
		dy_1 = dy_1[m1:m2]
		dz_1 = dz_1[n1:n2]
		dx_tilde = dx_tilde[l1:l2]
		dy_tilde = dy_tilde[m1:m2]
		dz_tilde = dz_tilde[n1:n2]
		ldegenerated = [ l1, m1, n1 ] eq [ l2, m2, n2 ]
		name += "trimmed_"
	endif

	if (not keyword_set (quiet)) then begin
		print, ' t = ', t
		print, ''
	endif

	time = t
	tags.time = t
	name += strtrim (xs, 2)+"_"+strtrim (xe, 2)+"_"+strtrim (ys, 2)+"_"+strtrim (ye, 2)+"_"+strtrim (zs, 2)+"_"+strtrim (ze, 2)
	sub_grid = create_struct (name=name, $
		['t', 'x', 'y', 'z', 'dx', 'dy', 'dz', 'Lx', 'Ly', 'Lz', 'dx_1', 'dy_1', 'dz_1', 'dx_tilde', 'dy_tilde', 'dz_tilde', 'lequidist', 'lperi', 'ldegenerated', 'x_off', 'y_off', 'z_off'], $
		t, x, y, z, dx, dy, dz, Lx, Ly, Lz, dx_1, dy_1, dz_1, dx_tilde, dy_tilde, dz_tilde, lequidist, lperi, ldegenerated, xns-nghostx, yns-nghosty, zns-nghostz)

	if (addghosts) then begin
		xs = xns - nghostx
		xe = xne - nghostx
		ys = yns - nghosty
		ye = yne - nghosty
		zs = zns - nghostz
		ze = zne - nghostz
	end

END


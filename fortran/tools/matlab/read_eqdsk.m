% reads the content of an EQDSK file (all in MKSA units)
%	function [out]=read_eqdsk(flnm,flpth)
% Inputs:
%		flnm	file name
%		flpth file path, if no path is entered, the default path is used
% This routine is not written at all in a general way, therefore it only works for a specific output file structure

function [out]=read_eqdsk(flnm,flpth)

 
% default path
if ~exist('flpth')||isempty(flpth)
%	flpth='/home/space/phsgbe/runs/chease/output/';
	flpth='./';
end

% reads file
if unix(['test -e ' flpth flnm])==0
	fid = fopen([flpth flnm], 'r');
else	
	error(['The file ' flpth flnm ' does not exist' ])
end
frewind(fid);

% comments and date
out.comments=fscanf(fid,'%c',48);

% switch
out.switch=fscanf(fid,'%i5');
%if out.switch~=3,
%	disp(['Warning, i3=' num2str(out.switch) ', unknown file structure...'])
%	return
%end

% box size
out.nrbox=fscanf(fid,'%5i',1);		% number of R points
out.nzbox=fscanf(fid,'%5i',1);	% number of Z points

% first line
out.rboxlength=fscanf(fid,'%e',1);	% bow width 
out.zboxlength=fscanf(fid,'%e',1);	% box height
out.R0EXP=fscanf(fid,'%e',1);	% Normalizing radius in CHEASE
out.rboxleft=fscanf(fid,'%e',1);	% left position of the R grid
out.zboxshift=fscanf(fid,'%e',1); % vertical displacement of the centre of the Z grid

% second line
out.Raxis=fscanf(fid,'%e',1);	% Raxis
out.Zaxis=fscanf(fid,'%e',1);	% Zaxis
out.psiaxis=fscanf(fid,'%e',1);	% psi_axis-psi_edge
out.psiedge=fscanf(fid,'%e',1);	% psi_edge-psi_edge (=0)
out.B0EXP=fscanf(fid,'%e',1);	% Normalizing magnetic field in CHEASE

% third line
out.current=fscanf(fid,'%e',1);	% total plasma current
fscanf(fid,'%e',4);% nothing or already stored

% fourth line
fscanf(fid,'%e',5);% nothing or already stored

% T(=RBt)
out.T=fscanf(fid,'%e',out.nrbox);

% p (pressure)
out.p=fscanf(fid,'%e',out.nrbox);

% TT'
out.TTprime=fscanf(fid,'%e',out.nrbox);

% p'
out.pprime=fscanf(fid,'%e',out.nrbox);

% psi
out.psi=fscanf(fid,'%e',out.nrbox*out.nzbox);
out.psi=reshape(out.psi,out.nrbox,out.nzbox)';

% q
out.q=fscanf(fid,'%e',out.nrbox); % (safety factor)

% nb points for the LCFS and limiter boundary
out.nLCFS=fscanf(fid,'%5i',1);
out.nlimits=fscanf(fid,'%5i',1);

% RZ LCFS description
tmp=fscanf(fid,'%e',out.nLCFS*2);
out.R_LCFS=tmp(1:2:end);
out.Z_LCFS=tmp(2:2:end);

% RZ limits
if out.nlimits>0
  tmp=fscanf(fid,'%e',out.nlimits*2);
  out.R_limits=tmp(1:2:end);
  out.Z_limits=tmp(2:2:end);
end


% RZ grid for psi
out.R_grid=linspace(out.rboxleft,out.rboxleft+out.rboxlength,out.nrbox);
out.Z_grid=out.zboxshift+linspace(-out.zboxlength/2,out.zboxlength/2,out.nzbox);


% psi grid for radial profiles
out.psi_grid=linspace(out.psiaxis,out.psiedge,out.nrbox);

% corresponding rhopsi
out.rhopsi=sqrt(abs(out.psi_grid-out.psiaxis)./abs(out.psiedge-out.psiaxis));

fclose(fid);

%% ########################################################################
%
%  
%  Version: 1.0
%
%  INPUT:
%
%  OUTPUT:
%
%  HISTORY:
%      v1.0: First runable state
%
%% ########################################################################
%% INITILIZATION
%--------------------------------------------------------------------------
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
format shortg;
%format compact;
set(0,'DefaultFigureWindowStyle','docked')

% Change the current folder to the folder of this m-file.
if(~isdeployed)
	cd(fileparts(which(mfilename)));
end
%%

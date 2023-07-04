% makes a simple plot of the directional spectra

try 
    load CR10_dspec_BDM.mat
catch 
    load CR10_dspec_BDM_supplied_output
    disp("This is using the supplied output")
    disp("run BDM_using_Matsuba_Code.m to generate your own")
end

figure(1); clf

timestep = 12;

pcolor(dspec.f, dspec.direction_nautical, dspec.S_f_theta_BDM(:,:,timestep)')

xlabel("Frequency [Hz]")
ylabel("Direction [from]")
c = colorbar;
    c.Label.String = "Energy Density [m^2/Hz/deg]";
dspec.dtime.Format = 'dd-MMM-uuuu HH:mm';
title(string(dspec.dtime(timestep)))


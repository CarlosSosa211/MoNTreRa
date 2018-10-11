tumDens = load('../OutputFiles/tumDens.dat');
itumDens = trapz(tumDens(:,1), tumDens(:,2))
itumDens = simpson(tumDens(:,1), tumDens(:,2))
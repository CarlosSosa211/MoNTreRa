clear all
close all
global nPar
global allTissues
global densTissues nonDensTissues;
global vascTissues nonVascTissues;
global b color nfig shape;
global fileNames outputNames;

densTissues = [1, 2, 5, 6, 8, 9, 11, 12, 19, 20, 21];
nonDensTissues = [3, 4, 7, 10, 13, 14, 15, 16, 17, 18];
vascTissues = [4, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20];
nonVascTissues = [1, 2, 3, 5, 6, 9, 15, 17, 19, 21];
allTissues = [densTissues, nonDensTissues];

% densTissues = [1, 2, 5, 6, 8, 9, 11, 19, 20];
% nonDensTissues = [3, 4, 7, 10, 14, 15, 16, 17];
% vascTissues = [4, 7, 8, 10, 11, 14, 16, 20];
% nonVascTissues = [1, 2, 3, 5, 6, 9, 15, 17, 19];
% allTissues = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 19, 20];

nDensTissues = length(densTissues);
nVascTissues = length(vascTissues);
nTissues = length(allTissues);

% path = '../../Carlos/Results/Morris100_39Par_Cluster';
% b = {'$tum$', '$T_{tum}$', '$N$', '$hea$', '$T_{heal}$', '$ang$'...
%     '$T_{end}$', '$D^{VEGF}$', '$V_{max}^{VEGF}$', '$K_M^{VEGF}$'...
%     '$\bar{v}$', '$v^{hyp}$', '$\alpha_{heal}$', '$\alpha/\beta_{heal}$'...
%     '$\alpha_{tumG1}$', '$\alpha/\beta_{tumG1}$', '$\alpha_{tumS}$'...
%     '$\alpha/\beta_{tumS}$', '$\alpha_{tumG2}$'...
%     '$\alpha/\beta_{tumG2}$', '$\alpha_{tumM}$', '$\alpha/\beta_{tumM}$'...
%     '$\alpha_{tumG0}$', '$\alpha/\beta_{tumG0}$', '$\alpha_{preEnd}$'...
%     '$\alpha/\beta_{preEnd}$', '$\alpha_{neoEnd}$'...
%     '$\alpha/\beta_{neoEnd}$', '$d$', '$d_{thres}$', '$T_{arrest}$'...
%     '$oxy$', '$pO_2^{nec}$', '$D^{O_2}$', '$V_{max}^{O_2}$'...
%     '$K_M^{O_2}$', '$pO_2^{preEnd}$', '$pO_2^{neoEnd}$', '$pO_2^{hyp}$'};
% bNum = ['tum (1), T_tum (2), N (3), hea (4), T_heal (5), ang (6),\n'...
%     'T_end (7), D^VEGF (8), V_max^VEGF (9), K_M^VEGF (10), barv (11), '...
%     'v^hyp (12),\nalpha_heal (13), alpha/beta_heal (14),'...
%     'alpha_tumG1 (15), alpha/beta_tumG1 (16),\nalpha_tumS (17), '...
%     'alpha/beta_tumS (18), alpha_tumG2 (19), '...
%     'alpha/beta_tumG2 (20),\nalpha_tumM (21), alpha/beta_tumM (22), '...
%     'alpha_tumG0 (23), alpha/beta_tumG0 (24),\nalpha_preEnd (25), '...
%     'alpha/beta_preEnd (26), alpha_neoEnd (27), '...
%     'alpha/beta_neoEnd (28),\nd (29), d_thres (30), T_arrest (31), '...
%     'oxy (32), pO_2^nec (33), D^O_2 (34),\nV_max^O_2 (35), '...
%     'K_M^O_2 (36), pO_2^preEnd (37), pO_2^neoEnd (38), pO_2^hyp (39)'];

path = '../../Carlos/Results/Morris100_38Par_2Gy_Cluster';
b = {'$N$', '$tum$', '$T_{tum}$', '$hea$', '$T_{heal}$', '$ang$'...
    '$T_{end}$', '$D^{VEGF}$', '$V_{max}^{VEGF}$', '$K_M^{VEGF}$'...
    '$\bar{v}$', '$v^{hyp}$', '$\alpha_{heal}$', '$\alpha/\beta_{heal}$'...
    '$\alpha_{tumG1}$', '$\alpha/\beta_{tumG1}$', '$\alpha_{tumS}$'...
    '$\alpha/\beta_{tumS}$', '$\alpha_{tumG2}$'...
    '$\alpha/\beta_{tumG2}$', '$\alpha_{tumM}$', '$\alpha/\beta_{tumM}$'...
    '$\alpha_{tumG0}$', '$\alpha/\beta_{tumG0}$', '$\alpha_{preEnd}$'...
    '$\alpha/\beta_{preEnd}$', '$\alpha_{neoEnd}$'...
    '$\alpha/\beta_{neoEnd}$', '$d_{thres}$', '$T_{arrest}$'...
    '$oxy$', '$pO_2^{nec}$', '$D^{O_2}$', '$V_{max}^{O_2}$'...
    '$K_M^{O_2}$', '$pO_2^{preEnd}$', '$pO_2^{neoEnd}$', '$pO_2^{hyp}$'};
bNum = ['N (1), tum (2), T_Tum (3), hea (4), T_heal (5), ang (6),\n'...
    'T_end (7), D^VEGF (8), V_max^VEGF (9), K_M^VEGF (10), barv (11), '...
    'v^hyp (12),\nalpha_heal (13), alpha/beta_heal (14), '...
    'alpha_tumG1 (15), alpha/beta_tumG1 (16),\nalpha_tumS (17), '...
    'alpha/beta_tumS (18), alpha_tumG2 (19), '...
    'alpha/beta_tumG2 (20),\nalpha_tumM (21), alpha/beta_tumM (22), '...
    'alpha_tumG0 (23), alpha/beta_tumG0 (24),\nalpha_preEnd (25), '...
    'alpha/beta_preEnd (26), alpha_neoEnd (27), '...
    'alpha/beta_neoEnd (28),\nd_thres (29), T_arrest (30), oxy (31), '...
    'pO_2^nec (32), D^O_2 (33), V_max^O_2 (34),\nK_M^O_2 (35), '...
    'pO_2^preEnd (36), pO_2^neoEnd (37), pO_2^hyp (38)'];

% path = uigetdir('../../Carlos/Results');
% path = '../../Carlos/Results/Morris100_34Par_Cluster';

nPar = length(b);
color = linspace(0, 1, nTissues);
shape = ['o', 's', 'v', 'd'];

nfig = 0;
quit = 0;

fileNames = {'EndTreatTumDens', '3MonTumDens'...
    'FinTumVol', 'IntTumDens', 'Killed50'...
    'Killed80', 'Killed90', 'Killed95', 'Killed99'...
    'Killed999', 'TimeTo95', 'TimeTo99', 'Rec'...
    'RecTumDens', 'RecTime'};
outputNames = {'Tumor density at the end of treat.'...
    'Tumor density 3 months after the beginning of treat.'...
    'Final tumor volume', 'Integral of tumor density'...
    '50% of tumor cells killed', '80% of tumor cells killed'...
    '90% of tumor cells killed', '95% of tumor cells killed'...
    '99% of tumor cells killed', '99.9% of tumor cells killed'...
    'Time to kill 95% of tumor cells'...
    'Time to kill 99% of tumor cells', 'Recurrence'...
    'Tumor density at recurrence', 'Recurrence time'};
nOut = 15;

while(~quit)
    selOut = input(['Select an output [endTreaTumDens (1), '...
        '3MonTumDens (2), finTumVol (3),\n intTumDens (4), '...
        '50%killed (5), 80%killed (6), 90%killed (7), 95%killed (8),\n '...
        '99%killed (9), 99.9%killed (10), timeTo95 (11), '...
        'timeTo99 (12), rec(13),\n recTumDens (14), recTime (15), '...
        'or all of them (-1)]: ']);

    if(selOut >= 1 && selOut <= nOut)
            par = input(['Select one parameter (', bNum...
                ')\nor all of them (0): ']);
            tissueSet = input(['Define a set [all (1), dense (2), '...
                'non-dense (3), vascularised (4) or non-vascularised '...
                '(5)] tissues: ']);
            if(par >= 1 && par <= nPar)
                plotParOutput(path, par, tissueSet, selOut)
            elseif(par == 0)
                for i = 1:nPar
                    plotParOutput(path, i, tissueSet, selOut)
                end
            end
                
    elseif(selOut == 0)
            quit = 1;
    end
end

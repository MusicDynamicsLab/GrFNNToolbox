
%% Make short cochlea network
nbm = networkMake(1, 'hopf', -0.1,       0, 0, 0, 0, 0.0000, ...
                     'log', 30, 5000, 320, 'channel', 1, ...
                     'display', 0, 'save', 10);
nbm.df = 1;
noc = networkMake(2, 'hopf',  0.00, -200000, 0, 0, 0, 0.0000, ...
                     'log', 30, 5000, 320, ...
                     'display', 0, 'save', 10);

%% Make CN and IC networks
ncn = networkMake(3, 'hopf', -1.0, 0, -100, 0, 0, 1.0, ...
                     'log', 10, 2500, 280, ...
                     'display', 0, 'save', 10);
nic = networkMake(4, 'hopf', -1.0, 0, -100, 0, 0, 1.0, ...
                     'log', 10, 2500, 280, ...
                     'display', 0, 'save', 10);

%% Add BM->OC connections
bm2oc = connectMake(nbm, noc, 'one', 1, 1, 0, 1);
noc   = connectAdd(nbm, noc, bm2oc, 'weight', 1, 'type', '1freq');

%% Somewhat fitted cochlear params, see GrFFNCochlea for latest parameter fits
nbm.a = -500 + i*2*pi.*nbm.f;
noc.b1 = -37000000;
noc.con{1}.w = 2000;

%% Add brainstem connections
oc2cn  = connectMake(noc, ncn, 'full',  1, 0, 1, 0);
cn2ic  = connectMake(ncn, nic, 'full',  1, 0, 1, 0);

ncn   = connectAdd(noc, ncn, oc2cn, 'weight', .5);
nic   = connectAdd(ncn, nic, cn2ic, 'weight', .1);

%% Run the network
M = modelMake(@zdot, @cdot, s, nbm, noc, ncn, nic);

tic;
M = odeRK4fs(M,s);
toc;

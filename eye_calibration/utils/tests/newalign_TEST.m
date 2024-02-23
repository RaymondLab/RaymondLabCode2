load('newalign_TEST.mat', 'mag1', 'mag2', 'vid', 'linearity_shifts');

vidPos = vid.pos_data_upsampled;
vidVel = vid.vel_data_upsampled;

mag1Pos = mag1.pos_data;
mag1Vel = mag1.vel_data;
mag2Pos = mag2.pos_data;
mag2Vel = mag2.vel_data;

a = mag2.pos_data;
b = vidPos;


a = a - mean(a);
b = b - mean(b);

shift2 = finddelay(a, b, 5000);

%a = butterworthfilter(a, 15, 1000, 9);
%b = butterworthfilter(b, 15, 1000, 9);

[cors, lags] = xcorr(b, a, 5000);
cors = cors / max(cors);
[maxval, maxid] = max(cors);

shift = lags(maxid);

% a= ; b= ; a=a-mean(a); b=b-mean(b); testShift=finddelay(a,b,5000);
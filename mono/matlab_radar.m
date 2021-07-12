% Set the seed for repeatability
rng(2021);

% System parameters
c = physconst('lightspeed');
fc = 1e9; % Hz
bandwidth = 1.5e6; % Hz
prf = 4e3; % Hz
updateRate = prf;
rangeRes = c/(2*bandwidth); % m
rangeAmb = c/(2*prf); % m
rangeLims = [0 2*rangeAmb];


% Scanning parameters
azLimits = [-20 20];
elLimits = [0 0];
fov = [4;8]; % deg
numFullScans = 20; % number of full scans to simulate

numScanPointsAz = floor(diff(azLimits)/fov(1)) + 1;
numScanPointsEl = floor(diff(elLimits)/fov(2)) + 1;
numScanPoints = numScanPointsAz*numScanPointsEl;
angRes = fov/4;

refTgtRange = 20e3;
refTgtRCS = 0;
detProb = 0.9;
faRate = 1e-6;

pitch = elLimits(1); % mount the antenna rotated upwards

% Create radar object
sensorIndex = 1; % a unique identifier is required
radar = radarDataGenerator(sensorIndex,'UpdateRate',updateRate,...
    'DetectionMode','Monostatic','ScanMode','Electronic','TargetReportFormat','Detections','DetectionCoordinates','scenario',...
    'HasElevation',true,'HasINS',true,'HasRangeRate',false,'HasRangeAmbiguities',true,'HasFalseAlarms',true,...
    'CenterFrequency',fc,'Bandwidth',bandwidth,...
    'RangeResolution',rangeRes,'AzimuthResolution',angRes(1),'ElevationResolution',angRes(2),...
    'ReferenceRange',refTgtRange,'ReferenceRCS',refTgtRCS,'DetectionProbability',detProb,'FalseAlarmRate',faRate,...
    'RangeLimits',rangeLims,'MaxUnambiguousRange',rangeAmb,'ElectronicAzimuthLimits',azLimits,'ElectronicElevationLimits',elLimits,...
    'FieldOfView',fov,'MountingAngles',[0 pitch 0]);

scanplt = helperScanPatternDisplay(radar);
scanplt.makeOverviewPlot;


numResolutionCells = diff(rangeLims)*prod(fov)/(rangeRes*prod(angRes));

% Create scenario
stopTime = numFullScans*numScanPoints/prf;
scene = radarScenario('StopTime',stopTime,'UpdateRate',0);

numFrames = numFullScans*numScanPoints; % total number of frames for detection
expNumFAs = faRate*numFrames*numResolutionCells % expected number of FAs


% Create platforms
rdrPlat = platform(scene,'Position',[0 0 12]); % place 12 m above the origin
tgtPlat(1) = platform(scene,'Trajectory',kinematicTrajectory('Position',[38e3 6e3 10e3],'Velocity',[-380 -50 0]));
tgtPlat(2) = platform(scene,'Trajectory',kinematicTrajectory('Position',[25e3 -3e3 1e3],'Velocity',[-280 10 -10]));


% Mount radar to platform
rdrPlat.Sensors = radar;

% Set target platform RCS profiles
tgtRCS = [20, 4]; % dBsm
tgtPlat(1).Signatures = rcsSignature('Pattern',tgtRCS(1));
tgtPlat(2).Signatures = rcsSignature('Pattern',tgtRCS(2));

time = (0:numFrames-1).'/updateRate;
tgtPosTruth(:,:,1) = tgtPlat(1).Position + time*tgtPlat(1).Trajectory.Velocity;
tgtPosTruth(:,:,2) = tgtPlat(2).Position + time*tgtPlat(2).Trajectory.Velocity;
tgtRangeTruth(:,1) = sqrt(sum((tgtPosTruth(:,:,1)-rdrPlat.Position).^2,2));
tgtRangeTruth(:,2) = sqrt(sum((tgtPosTruth(:,:,2)-rdrPlat.Position).^2,2));





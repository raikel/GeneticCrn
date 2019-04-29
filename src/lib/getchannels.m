function mC = getchannels(options)

% Number of points
n = options.linksQ;

switch options.geometry
  case 'square'
    % Scenary Dimensions
    xLim = options.xLim;
    yLim = options.yLim;
    % Transmitters positions
    xPosTx = ( xLim(2) - xLim(1) )*rand(1, n) + xLim(1);
    yPosTx = ( yLim(2) - yLim(1) )*rand(1, n) + yLim(1);
    posTx = [xPosTx; yPosTx];    
  case 'cell'
    % Scenary Dimensions
    radius = options.radius;    
    % Transmitters positions
    rhoPosTx = radius*sqrt( rand(1, n) );
    thetaPosTx = 2*pi*rand(1, n);
    % Back to cartesian coordinate system
    xPosTx = rhoPosTx.*cos(thetaPosTx);
    yPosTx = rhoPosTx.*sin(thetaPosTx);
    posTx = [xPosTx; yPosTx];   
end

% Receivers positions

% Optional fix transmitters positions
if ~isempty(options.fixTxPos)
  idxFixTxPos = ~isnan(options.fixTxPos);
  posTx(idxFixTxPos) = options.fixTxPos(idxFixTxPos);
end
% Transmitting range
txRange = options.txRange;
% If transmitting range is not specified
if isempty(txRange)
  % Set the transmitting range to be the half of distance to the most
  % closest transmitter
  dist = distmatrix(posTx, posTx);
  dist(dist == 0) = Inf;
  txRange = min(dist, [], 1)/2;          
end
% The angle of vector from transmitter to receiver is obtained as a
% random uniformly distributed number between 0 and 2*PI
rxAngle = 2*pi*rand(1, n);    
% Set the receivers positions
posRx = posTx + [txRange.*cos(rxAngle); txRange.*sin(rxAngle)];
    
% Plot net if plot option is set
if (options.plotNet == true)
  % Create a new figure
  figure;
  % Draw deployment area
  switch options.geometry
    case 'square'
    case 'cell'    
      t = linspace(eps, 2*pi, 32);
      fill(radius*sin(t), radius*cos(t), 'w')
      span = radius + max(txRange);
      axis([-span span -span span])
      axis square
      axis off  
  end
  % Hold cuurent plot
  hold on
  % Plot links
  plot([posTx(1, :); posRx(1, :)], [posTx(2, :); posRx(2, :)], 'b')
  % Plot transmitters as green squares and recievers as red circles 
  plot(posTx(1, :), posTx(2, :), 's', 'MarkerFaceColor', 'g');
  plot(posRx(1, :), posRx(2, :), 'o', 'MarkerFaceColor', 'r');
    
  for k = 1:n
      text(posTx(1, k), posTx(2, k), ['   ', num2str(k)] ,'FontSize', 6)
  end
end

dist = getcrossdist(posTx, posRx); % receptores por filas y tx por columnas

mC = linkmodel(dist);

%==========================================================================
function mC = linkmodel(dist)

alpha = 3;

mC = dist.^(-alpha);
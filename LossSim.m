classdef LossSim < handle
    
    properties
        DegradeValues = .1:.1:.9;
        AOI = [0.33, 0.33, 0.66, 0.66]
    end
    
    properties (SetAccess = private)
        SimData
        SimIDs
        AOIScores_Original
        AOIScores_Sim
    end
    
    properties (Dependent)
        TopGoodN
        TopGoodProp
        GapSampleSpan        
    end
    
    properties (Dependent, SetAccess = private)
        NumSubjects
        Subjects
        LossTable
    end
    
    properties (Access = private)
        ids
        x
        y
        t
        gaps
        numGaps
        propMissing
        topGoodProp = .1
        gapSampleSpan = .05
        simMissing
    end
    
    methods
        
        function obj = LossSim
            obj.x = {};
            obj.y = {};
            obj.t = {};
            obj.gaps = {};
            obj.numGaps = [];
            obj.propMissing = [];
            obj.TopGoodProp = .1;
        end
        
        function AddTrial(obj, id, x, y, t)
            % check input args
            if isnumeric(id)
                id = num2str(id);
            elseif ~ischar(id)
                error('ID must be numeric of char.')
            end
            if ~isnumeric(x) || ~isvector(x) || ~isnumeric(y) ||...
                    ~isvector(y) || ~isnumeric(t) || ~isvector(t)
                error('x, y, and t must be numeric vectors.')
            end
            if ~isequal(size(x), size(y), size(t))
                error('x, y, and t must be of the same size.')
            end
            % add
            cursor = length(obj.x) + 1;
            obj.x{cursor} = x;
            obj.y{cursor} = y;
            obj.t{cursor} = t;
            obj.ids{cursor} = id;
            % calculate 
            [obj.gaps{cursor}, obj.numGaps(cursor),...
                obj.propMissing(cursor)] = etMeasureGaps(x, y);
        end
        
        function [data, gaps] = SubjectData(obj, id)
            % find id
            found = find(strcmpi(obj.ids, id));
            if isempty(found)
                data = [];
            end
            numTrials = length(found);
            data = cell(1, numTrials);
            for tr = 1:numTrials
                data{tr} = [...
                            obj.t{found(tr)},...
                            obj.x{found(tr)},...
                            obj.y{found(tr)}];
                gaps{tr} = obj.gaps{found(tr)};
            end
        end
        
        function Simulate(obj)
            wb = waitbar(0, 'Simulating');
            % get loss table
            lt = obj.LossTable;
            % loop through degrade values, get gaps from relevant datasest
            % matching that value (e.g. datasets with lost data between
            % 0.05 and 0.15 for value 0.1 with span 0.05). Fit a PDF to the
            % distribution of gaps, and store it for later use
            numSim = length(obj.DegradeValues);
            simVal = false(1, numSim);
            gapPdf = cell(numSim, 1);
            for v = 1:numSim
                wb = waitbar((v / numSim) * .1, wb, 'Fitting gap distributions');
                % get current value, find edges of span
                val = obj.DegradeValues(v);
                v1 = val - obj.GapSampleSpan;
                v2 = val + obj.GapSampleSpan;
                % find datasets that match these values
                idx = lt{:, 3} >= v1 & lt{:, 3} < v2;
                % only continue if there is at least one subject worth of
                % data
                if sum(idx) > 1
                    simVal(v) = true;
                    % get gaps from these datasets
                    [~, gaps] = cellfun(@(x) obj.SubjectData(x), lt{idx, 1},...
                        'uniform', false);
                    % collapse across trials
                    gaps = cellfun(@(x) vertcat(x{:}), gaps, 'uniform', false);
                    % collapse across subjects
                    gaps = vertcat(gaps{:});
                    % fit pdf
                    gapPdf{v} = fitdist(gaps, 'inverse gaussian');
                end
            end
            
            %  adjust number of sims to possible data. If a sim has been
            %  specified with 90% missing data, we may not have a dataset
            %  that with that high an amount missing. In this case we can't
            %  run this sim, so we exclude that value. 
            if ~all(simVal)
                warning('Not all degradation values could be used due to insufficient lost data. These values have been removed.')
                obj.DegradeValues(~simVal) = [];
                gapPdf(~simVal) = [];
                numSim = length(obj.DegradeValues);
            end
            
            % get subject data for top N good datasets
            ids = lt{1:obj.TopGoodN, 1};
            obj.SimIDs = ids;
            data = cellfun(@(x) obj.SubjectData(x), ids, 'uniform', false);
            data = cellfun(@(x) vertcat(x{:}), data, 'uniform', false);
            
            % for each sim and each subject, make a vector of simulated
            % missing values
            obj.simMissing = cell(obj.TopGoodN, numSim);
            obj.SimData = cell(obj.TopGoodN, numSim);
            for s = 1:obj.TopGoodN
                wb = waitbar(0.1 + ((s / obj.TopGoodN) * .9), wb,...
                    'Simulating data loss');
                for v = 1:numSim
                    missing = all(isnan(data{s}(:, 2:3)), 2);
                    numSamps = length(missing);
                    % get target prop of missing data
                    target = obj.DegradeValues(v);
                    % get actual missing data
                    actual = sum(missing) / numSamps;
                    % interatively remove data until target is reached
                    while actual < target
                        % sample a gap duration, if gap duration is longer
                        % than num samples, resample
                        gapdur = round(random(gapPdf{v}));
                        while gapdur > numSamps
                            gapdur = round(random(gapPdf{v}));
                        end
                        % choose a random onset
                        onset = randi(numSamps - gapdur);
                        % insert random data
                        missing(onset:onset + gapdur) = true;
                        % recalc loss
                        actual = sum(missing) / numSamps;
                    end
                    % store missing values
                    obj.simMissing{s, v} = missing;
                    % copy data, inset missing values
                    obj.SimData(s, v) = data(s);
                    obj.SimData{s, v}(missing, 2:3) = nan;
                end
            end
            obj.ProcessAOIs;
            delete(wb)
        end
        
        function ProcessAOIs(obj)
            if obj.NumSubjects == 0
                error('No subjects.')
            end
            % preallocate
            obj.AOIScores_Original = nan(obj.TopGoodN, 1);
            obj.AOIScores_Sim =...
                nan(obj.TopGoodN, length(obj.DegradeValues));
            % process each subject in turn
            for s = 1:obj.TopGoodN
                % original data
                data = obj.SubjectData(obj.Subjects{s});
                data = vertcat(data{:});
                inaoi_orig = obj.inaoi(data(:, 2), data(:, 3));
                obj.AOIScores_Original(s) =...
                    sum(inaoi_orig) / sum(all(~isnan(data(:, 2:3)), 2));
                % simulated data
                if ~isempty(obj.SimData)
                    inaoi_sim = cellfun(@(x) obj.inaoi(x(:, 2), x(:, 3)),...
                        obj.SimData(s, :), 'uniform', false);
                    obj.AOIScores_Sim(s, :) = cellfun(@(x, y)...
                        sum(x) / sum(all(~isnan(y(:, 2:3)), 2)),...
                        inaoi_sim, obj.SimData(s, :));
                else
                    obj.AOIScores_Sim = [];
                end
            end
        end
        
        function PlotScores(obj)
            if obj.NumSubjects == 0
                error('No subjects.')
            end
            % start drawing
            figure
            subplot(4, 1, 1:3)
            hold on
            numSim = length(obj.DegradeValues);
            x = obj.DegradeValues; 
            dy = nan(obj.TopGoodN, numSim - 1);
            dx = diff(x);
            if ~isempty(obj.AOIScores_Sim)
                % calculate and draw the mean, sd and sem
                mu = mean(obj.AOIScores_Sim, 1);
                sd = std(obj.AOIScores_Sim, 1);
                sem = sd / sqrt(obj.TopGoodN);
                plot(x, mu, 'Marker', 'x',...
                    'MarkerSize', 10, 'MarkerEdgeColor', 'k',...
                    'Color', 'k');
                plot(x, mu + (2 * sem), '--', 'Color', 'k');      
                plot(x, mu - (2 * sem), '--', 'Color', 'k');   
                plot(x, mu + (2 * sd), '-.', 'Color', 'r');      
                plot(x, mu - (2 * sd), '-.', 'Color', 'r');    
                % axis details
                xlabel('Simulated prop. data loss')
                ylabel('Proportion looking time in AOI')
                xl = xlim;
%                 legend(obj.SimIDs, 'location', 'BestOutside')
                % draw each subject 
                for s = 1:obj.TopGoodN
                    y = obj.AOIScores_Sim(s, :);
                    dy(s, :) = diff(y);
                    pl = plot(x, y, 'Marker', 'o');
                    pl.Color(4) = .2;
                    scatter(x, y)
                end
            end
            ylim([0, 1])
            legend('Mean', '+2SEM', '-2SEM', '+2SD', '-2SD');
            % calculate change in slope between values
            subplot(4, 1, 4)
            area(x, [nan, mean(abs(dy), 1)]);
            hold on
            plot(x, sd - min(sd), 'y', 'linewidth', 3)
            xlim(xl);
            ylabel('Mean change in prop. looking time')
            legend('Prop looking change', 'SD change', 'location', 'best')
        end
        
        % get / set methods
        function val = get.LossTable(obj)
            % get list of subjects
            subs = obj.Subjects;
            numSubs = length(subs);
            numGaps = nan(numSubs, 1);
            propMissing = nan(numSubs, 1);
            % loop through and get all trials
            for s = 1:length(subs)
                found = strcmpi(obj.ids, subs{s});
                numGaps(s) = sum(obj.numGaps(found));
                propMissing(s) = mean(obj.propMissing(found));
            end
            val = cell2table([subs', num2cell(numGaps), num2cell(propMissing)],...
                'variablenames', {'ID', 'NumGaps', 'PropMissing'});
            % rank
            [~, val.Rank] =...
                ismember(val.PropMissing, unique(val.PropMissing));
            % sort by rank
            val = sortrows(val, 'Rank');
        end
        
        function val = get.Subjects(obj)
            val = unique(obj.ids);
        end
        
        function val = get.NumSubjects(obj)
            val = length(obj.Subjects);
        end
        
        function val = get.TopGoodN(obj)
            if obj.NumSubjects == 0
                val = [];
            else
                val = ceil(obj.NumSubjects * obj.topGoodProp);
            end
        end
        
        function set.TopGoodN(obj, val)
            if obj.NumSubjects == 0
                error('Cannot set TopGoodN when no data has been added.')
            end
            % check args
            if ~isnumeric(val) || ~isscalar(val) || val < 1
                error('Value must be a positive scalar.')
            end
            if val > obj.NumSubjects
                error('Value exceeds number of subjects.')
            end
            % calculcate percentage
            obj.topGoodProp = val / obj.NumSubjects;
        end
        
        function val = get.TopGoodProp(obj)
            val = obj.topGoodProp;
        end
        
        function set.TopGoodProp(obj, val)
            % check args
            if ~isnumeric(val) || ~isscalar(val) || val <= 0 || val > 1
                error('Value must be a proportion between 0 and 1.')
            end
            obj.topGoodProp = val;
        end
        
        function set.GapSampleSpan(obj, val)
            % check that span does not exceed degrade values
            if obj.DegradeValues(1) - val < 0 || val + obj.DegradeValues(end) > 1
                error('GapSampleSpan must not exceed upper or lower bounds of DegradeValues.')
            end
            obj.gapSampleSpan = val;
        end
        
        function val = get.GapSampleSpan(obj)
            val = obj.gapSampleSpan;
        end
        
        function set.AOI(obj, val)
            obj.AOI = val;
            obj.ProcessAOIs;
        end

    end
    
    methods (Access = private)
        
        function in = inaoi(obj, x, y)
            in = x >= obj.AOI(1) & x <= obj.AOI(3) &...
                y >= obj.AOI(2) & y<= obj.AOI(4);
            in = double(in);
        end
        
    end
end
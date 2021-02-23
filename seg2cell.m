function ls = seg2LossSim

    % get unique IDs
    [id_u, ~, id_s] = unique(seg.ids);
    id_num = length(id_u);
    
    % loop through IDs and put in trial data
    for i = 1:id_num
        
        % get seg entries for this subject
        segIdx = find(id_s == i);
        
        % get number of trials, and max number of samples for these trials
        numTrials = length(segIdx);
        
        % loop through and place data in c
        col = 1;
        for t = 1:numTrials(i)
            c{i}(1:numRows{i}(t), col) = seg.x{segIdx(t)};
            c{i}(1:numRows{i}(t), col + 1) = seg.y{segIdx(t)};
            col = col + 2;
        end
        
    end

end
function [combinations] = get_pairwise_combs(CELL_TYPEs,NUMMEL_BRs)

% Get field names and values
cellTypes = fieldnames(NUMMEL_BRs);
cellValues = struct2cell(NUMMEL_BRs);

% Initialize empty cell array to store combinations
combinations = cell(0, 4);
% Generate pairwise combinations
for i = 1:numel(cellTypes)
    if size(cellValues{i},2)>1
        pairs = nchoosek(cellValues{i},2);
        combinations = [combinations; repmat({cellTypes{i}}, size(pairs, 1), 1), num2cell(pairs(:, 1)), ...
            repmat({cellTypes{i}}, size(pairs, 1), 1), num2cell(pairs(:, 2))];
    else
        pairs = 1;
    end

    for j = (i+1):numel(cellTypes)
        [X, Y] = meshgrid(cellValues{i},cellValues{j});
        pairs = [X(:) Y(:)];

        combinations = [combinations; repmat({cellTypes{i}}, size(pairs, 1), 1), num2cell(pairs(:, 1)), ...
            repmat({cellTypes{j}}, size(pairs, 1), 1), num2cell(pairs(:, 2))];
    end
end


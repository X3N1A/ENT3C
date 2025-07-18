function [comparisons] = get_pairwise_combs(SAMPLES)

SAMPLE_IDX = 1:length(SAMPLES);
if length(SAMPLE_IDX)>1
    comparisons = nchoosek(SAMPLE_IDX,2);

    comparisons = arrayfun(@(x) SAMPLES(x), comparisons, 'UniformOutput', false);
else
    comparisons=SAMPLES;

end

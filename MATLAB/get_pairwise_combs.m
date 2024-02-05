function [comparisons] = get_pairwise_combs(SAMPLES)

SAMPLE_IDX = 1:length(SAMPLES);

comparisons = nchoosek(SAMPLE_IDX,2);

comparisons = arrayfun(@(x) SAMPLES(x), comparisons, 'UniformOutput', false);
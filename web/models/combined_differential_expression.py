"""
make an endpoint that can looks for common DEGs across metadata (e.g. organs, diseases. or cell types)

at first iteration, user supply disease of interest and a list of tissues

use disease and tissue as keywords to find datasets that satistfdy ther condition (containing both the disease and the tissue)
can use similar logic as existing code

the for each tissue:
    compute top DE genes (take inspiration from existing code for getting the top DE genes for a given disease and cell-type)
    more explanation on this:
        see the input for compute_diff_expression? you can supply a filter to the function at the end
        the filter is expected to be:
            {
                "disease": "xxx",
                "celltype": "xxx",
                "sex": "xxx",
                "tissue": "xxx", 
                ...
            }
        
        In this case we are going to do something similar, since we are only looking at one tissue in one disease
        we expect the filter to look like this"
            {
                "disease": "xxx",
                "tissue": "xxx",
            }
        
        when filtering the dataset, we first filter on the disease.
        then we filter on tissue.
        
        The rest should be similar to compute_diff_expression
        
        What should be returned from this function you ask?
        Based on the existing return structure from compute_diff_expression, 
        you would probably want to add things like which tissue. 
        We are still comparing between disease and normal, but just need to specify in what tissue. 
        
        At first, this will be averaged across all cell types (because the user does not need to supply a cell type filer).
        So we don't need to specify the cell_type. Only add cell type later when we allow the user to specify in that cell types (a list of cell types). 
        
        
        
        
    
we now have a list of top de genes for each disease-tissue pair
return twin intersection, triplet intersection, etc 

twin intersection = for any TWO, common DEs = [...]
triplet intersection = for any THREE, common DEs = [...]
(see the paper figure 5A, look at how it considered 1 tissue at a time, 2 tissues at a time, 3 tissues at a time)
(you will expect to see less genes when more tissues are considered together)


"""
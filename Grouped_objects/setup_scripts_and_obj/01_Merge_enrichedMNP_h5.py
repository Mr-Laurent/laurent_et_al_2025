#!/usr/bin/env python3
import anndata as ad
import os

def main():
    os.chdir('./Metacells_EnrichMNP/')

    adn2t=ad.AnnData.transpose(ad.read_h5ad('./umitabcr5_n2.h5ad'))
    adn4t=ad.AnnData.transpose(ad.read_h5ad('./umitabcr5_n4.h5ad'))
    adn8t=ad.AnnData.transpose(ad.read_h5ad('./umitabcr5_n8.h5ad'))
    adn18t=ad.AnnData.transpose(ad.read_h5ad('./umitabcr5_n18.h5ad'))
    adn34t=ad.AnnData.transpose(ad.read_h5ad('./umitabcr5_n34.h5ad'))
    adn48t=ad.AnnData.transpose(ad.read_h5ad('./umitabcr5_n48.h5ad'))
    adn50t=ad.AnnData.transpose(ad.read_h5ad('./umitabcr5_n50.h5ad'))
    
    adall = adn2t.concatenate(
        adn4t,adn8t,adn18t,adn34t,adn48t,adn50t,
        join='inner',
        batch_categories=['adn2t','adn4t','adn8t','adn18t','adn34t','adn48t','adn50t'],
        uns_merge="first",
        index_unique='-'
    )

    ad.AnnData.write(ad16nov,filename="ad16nov22_reAug24.h5ad")

if __name__ == "__main__":
    main()
import numpy as np
from sklearn.manifold import TSNE
import sys

chr=sys.argv[1]
layer_output = np.load('/users/mtaranov/dnaseEmbed_blood/data/last_layer_embedding_chr'+chr+'.npy')
print "Computing TSNE comps..."
model = TSNE(n_components=2, random_state=0)
#np.set_printoptions(suppress=True)
first2comp = model.fit_transform(layer_output) 

print first2comp.shape

np.save('/users/mtaranov/dnaseEmbed_blood/data/'+'tsne_chr'+chr+'.npy', first2comp)


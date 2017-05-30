# Absorption Image Processing

The functions included in this git repository are designed to convert pictures taken in absorption imaging into useful images of atoms.  However as general image processing functions they may be helpful for other purposes.  Abstractly, it's general purpose is to reconstruct areas of an image using part of that image and a database of similar images.  Feel free to copy/edit/redistribute this code.  If you use do use it, I'd be happy to hear about how well, or how badly, it's worked for you.  Shoot me an email at zakven{at}mit{dot}edu.  While developing the code I kept in mind that it may be useful for others, so I tried to document it well and tried to not make it too lab-specific (e.g. by assuming specific naming conventions, etc.).  That being said, I did still make it somewhat specific as I was somewhat lazy not completely convinced that anyone else would use it.  If enough people express interest in the code, I will go back and try to make it more general and useful to others.

## Brief Description of Absorption Imaging

Absorption imaging is widely used in atomic physics.  Clouds of atoms are prepared, then a laser beam is shot through the atom cloud to a camera.  The projected 2D density distribution of the atoms can then be calculated by looking at the shadow of the semi-transparent atomic cloud.  Typically this process is done by taking two or three images.  First an image of the beam with the atomic cloud shadow is taken, which we'll call the raw image (in our lab we save these with the suffix "\_raw.ascii").  Next, the atom cloud is allowed to fly away and an image of the beam without any shadow is taken, which we'll call the background image (we save these with the suffix "\_back.ascii").  Optionally a third image is taken with the imaging laser completely off, which we'll call the noise image (not surprisingly we save these with the suffix "\_noise.ascii").  The naming convention for the last two is a little confusing, I think instead of "background" and "noise" they'd be better named "beam" and "background" respectively.  However, to remain consistent with the previous convention I'll use the old naming scheme here instead.  If the noise image is taken, the first step is to subtract it from both the raw and background images.  Once that is done, the optical depth of the cloud along the imaging beam can be calculated as OD=log( abs( back./raw ) ) where "./" means element-wise division and the absolute value is taken to avoid taking the log of a negative number.  This technique is implemented in get\_OD\_simple().  Typically many such optical depth measurements are taken and averaged to reduce statistical noise.

In practice this simple algorithm is complicated by the fact that interference patterns typically appear in the imaging beam.  Those fringe patterns often have systematic effects (e.g. due to the timing of shutters) and long term drifts due to temperature fluctuations.  Those systematic fringe patterns remain in the image even after averaging many shots.  A quick and dirty way to get rid of them is to do many runs of the experimental sequence but without making an atom cloud (e.g. block the MOT laser so no atoms are captured, but keep running the sequence so that shutter movement still systematically affects the images), and then process and average those images.  This will give the residual systematic intereference pattern, which can then be subtracted from the images with atoms.  However, this quick-and-dirty technique requires wasting time by running the experiment without any atoms.  A more clever background subtraction algorithm that doesn't require extra runs of the exeriment without atoms would save time.  Two such algorithms are included here and discussed below.

Experimentally we've seen that the two techniques give very similar images to each other, and those images are also similar to the results of the quick-and-dirty method discussed here.  The eigenfaces technique turns out to be more computationally efficient than the svd method, so it is the recommended algorithm.

## Quick Start

The quickest way to learn to use this code is probably by example.  Here is a typical usage.

We have a folder named "20170405" full of data from absorption imaging.  All the background image files end with "\_back.ascii" and we'd like to use all of them to form a basis.  The atoms in all of the images are contained in rows 40 to 60 and columns 50 to 80 of the image; the rest of the image will be assumed to be free of atoms by the algorithm.  We use the code provided here to make a basis for the subspace of background images, use it to reconstruct and subtract the background from a raw image "Cool100d100d80PGCZ4.4\_1\_raw.ascii".  We then plot the results.

```matlab
>> %Get the image we'd like to analyze
>> filename = fullfile('20170405','Cool100d100d80PGCZ4.4_1_raw.ascii');
>> image_in = load_image(filename);
>> 
>> %Select a background region
>> row_min=40; row_max=60; col_min=50; col_max=80;
>> back_region = make_back_region(image_in,row_min,row_max,col_min,col_max);
>> 
>> %Make a basis
>> max_vectors = 20; %20 is typically a good number for this
>> ls_pattern = fullfile('20170405','\*_back.ascii');
>> file_list = get_file_list(ls_pattern);
>> [basis_eig, mean_back] = make_basis_eig(file_list,back_region,max_vectors);
>> 
>> %Use the basis to get the atomic cloud's optical depth
>> OD_eig = get_OD_eig(image_in,basis_eig,mean_back,back_region);
>> 
>> %Plot the results
>> plot_image(image_in,['Original ',filename]);
>> plot_image(OD_eig,'Eigenfaces Optical Depth');
```

To understand what is going on behind the scenes in this example, read the subsequent sections of the README and/or the docstrings of the functions themselves.

Note that it's not necessary to contruct a new basis each time an image is analyzed.  Many images can be analyzed in a row with the same basis.  Generally a new basis only needs to be generated when the imaging beam changes significantly from how it was in the input images to make\_basis\_eig().  How often the basis needs to be re-generated will depend on the lab and how widely varying the input images to make\_basis\_eig() are.

## Singular Value Decomposition Method

The general idea used by the singular value decomposition (svd) method here is somewhat common in atomic physics, and was mentioned to us by Cheng Chin from the University of Chicago.  At it's heart the concept is based on undergrad-level linear algebra.  Consider a vector space of m-by-n images.  Addition of two images is defined by simply adding the photon counts of each corresponding pixel together.  The inner product of two images is defined by multiplying corresponding pixels together, then summing all of them.  This is equivalent to taking each column of the image and stacking them on top of each other to form a column vector, then using the usual definition of vector addition and inner product.  In fact, the code does exactly that to take advantage of matlab's buitl-in numerical algorithms.

Now that we have a vector space defined, consider the series background images (images of just the imaging beam with no atoms).  Due to temperature drifts, vibrations, etc. they are not all exactly the same vector in the vector space.  In reality, there is a subspace of possible background images.  Because the images are all very similar, it is a relatively small subspace compared to the full vector space, which has dimension mn.  The idea is to create a basis for the subspace of possible background images, then use it to reconstruct the background of a raw image.  That reconstructed background can then be used as the reference background image and the formula OD=log( abs( back./raw ) ) can then be applied.  The technique to make that basis and then use it to reconstruct the background is discussed in this section.

The most straight-forward way to get a basis is to take a handful of background images and orthogonalize them using Graham-Schmidt, or a more numerically stable equivalent of it.  In our experience, exactly which background images you use is very important.  We found that the best results were only obtained when we used some background images taken before the raw image and some taken after (Although this is likely very lab-depenent and your mileage may vary).  This is most likely because subsequent images are very similar, they likely only differ due to small vibrational movements of the optics, while images taken a long time apart are more different, as temperature drifts can move the beam much larger distances.  It would be nice if the algorithm would just take a large batch (we typically use ~1000) of images and made a good basis that had only a handful of vectors that still varied enough to span the subspace well.  To do this, the code puts all the background images into an array (each column is an image) and then takes the singular value decomposition of that matrix.  Then the vectors with the largest singular values are kept and used as the basis, and the rest are discarded.

Now that the basis is formed, it can be used to reconstruct the background of a raw image.  Exactly how to account for the region of the image with the atomic cloud is nontrivial and is discussed in its own section below.  Here we just give a qualitative overview of the idea.  The region of the raw image without atoms is projected onto each vector of the basis.  We then recontruct a beam profile for the whole image by summing up the basis images weighted by their projections.  This is a typical projection into a subspace from linear algebra, except that we ignore the region with atoms.  The reconstructed image is then used as the reference background.

## Eigenfaces Algorithm

While looking for papers about the algorithm mentioned above, I came across the paper "Reduction of interference fringes in absorption imaging of cold atom cloud using eigenface method" by Li et al.  The authors there used the facial recognition algorithm from "Eigenfaces for Recognition" by Turk and Pentland and applied it to absorption imaging.  The code included here implements the algorithm from Li et al. with some modifications to properly account for the region of the image with atoms.

Here I'll provide a brief description of the algorithm for reference.  For a complete explanation of the algorithm see the cited papers and the comments in the code.  The basic idea is to do principal component analysis as follows.  First we take the average off all the input background images, and then subtract the average from each  background image.  This gives us a bunch of images showing the variations from the mean background.  We take those images and turn them into vectors by stacking their columns, as described in the svd section.  Those vectors are then put next to each other to form a matrix we'll call A.  The column space of A then (approximately) spans the subspace of all possible variations from the mean background.

Now suppose we take another image.  We can then subtract the mean background from it as a first step and call the result Phi.  Next we project Phi onto the column space of A, which will give an image in the subspace of possible variations from the mean background.  In the end we'll use that projection (plus the mean background image) as our reconstructed background.  Doing that projection into the column space of A is equivalent to multiplying that vector by AA' where A' is the transpose of A.  To understand how that projection acts, it's helpful to express the Phi in terms of the eigenvectors of AA'.  If an eigenvector has a large eigenvalue, then it will contribute significantly to the projection.  If the eigenvalue is tiny, then acting on that eigenvector with AA' will give an image that is almost zero, and so it will not significantly contribute to the projection.  Therefore we can cleverly pick a small number of basis vectors that contribute significantly to this projection by choosing the eigenvectors with the largest eigenvalues.  The other eigenvectors can be neglected.  In this way we can take in many background images and distill them down to a handful of images that still can reconstruct backgrounds very well.

The only problem with what's been described so far is that the matrix AA' is very large.  If there are M images and one image is j-by-k, then the matrix A will be jk-by-M and AA' will be jk-by-jk.  That makes the computation of eigenvectors computationally intensive.  However, the eigenfaces algorithm uses a clever trick to find eigenvectors more easily.  The matrix A'A is M-by-M, which is typically much smaller than AA'.  Now suppose that v is an eigenvector of A'A, which we can find relatively easily due to the small size of A'A.  This implies that A'Av=lambda\*v where lambda is the eigenvalue.  Now if we left-multiply that expression by A we see that AA'(Av)=lambda\*(Av).  That implies that Av is an eigenvector of AA' with the same eigenvalue lambda.  Therefore we can get M eigenvectors of AA' by finding the eigenvectors of A'A and multiplying them by A.  Furthermore, we can ditch the ones with the smallest eigenvalues.  These eigenvectors (once re-shaped from column vectors to 2D images) are called the eigenfaces, due to the initial application of this algorithm to facial recognition.  Since AA' has jk eigenvectors while A'A only has M eigenvectors, we can't get all of the eigenvectors of AA' this way.  However, according to Li et al. this technique actually gives the M eigenvectors of AA' with the largest eigenvalues.  Admittedly I'm not sure how to show that, so it's on the to-do list.  TODO: prove eigenvectors are the ones with the largest eigenvalues

The eigenfaces algorithm gives very similar output to the svd algorithm.  So similar in fact that I've wondered if they're symbolically identical, i.e. they would give the exact same answer if evaluated exactly.  After some thought, I believe they would be if it weren't for the fact that the mean background is subtracted in the eigenfaces algorithm but not in the svd algorithm.  The similarities are particularly clear if you generate a basis with each algorithm with the same set of input images.  The first basis vector from the svd basis looks like the first vector from the eigenfaces basis, but with mean background image added to it.  Subsequent basis images look almost identical between the two (without even adding the mean background), up to a possible multiplication by negative one.  In practice, the eigenfaces algorithm runs much faster that the svd algorithm for large input sets.  For ~5000 images 121-by-121 pixels, it took my computer ~30 seconds to make a basis with svd and less than one second to make it with eigenfaces.  In practice it actually took another 15 seconds for each just to read the images from the hard drive.  To help remedy this, the load\_image() function caches hard drive reads so that future look ups of the same data are faster.


## Accounting For The Cloud Region

In Li et al. the region of the image with atoms was set to zero before doing the projection.  This works in practice if the atoms only take up a small portion of the image, however it does introdue systematic issues when reconstructing the background.  The easiest way to see that this is the case is as follows.  Contstruct an image that is the sum of the mean background and e where e is one of the eigenfaces.  If we then process that image as if it were a picture of atoms, we'd expect to find that the optical depth is zero everywhere since this image should project exactly onto the background image subspace.  However, if we set a region of the image to zero, then subtract the mean background, we'll get a vector that is different than the eigenface e.  Therefore, the projection onto e will not be exactly one and the reconstructed background will not be exactly right.

The way to correct for this is to contruct a basis that is orthonormal over the background region.  This method works for both the svd algorithm and the eigenfaces algorithm.  The technique is essentially the following.  Set the atom region of the background images to zero, then use those images to construct an orthonormal basis.  Next figure out how to write that orthonormal basis as a superposition of the input images (with the atom region still set to zero).  Now take that same superposition of images, except without setting the atom region to zero.  This gives a basis that is not orthonormal in the usual sense, but is orthonormal over the background region in the sense that if the atom region is set to zero, the basis is orthonormal.

Now we can use that basis appropriately by doing the following.  First, set the atom region of the input image to zero, then calculate its projection onto each basis vector.  Since the atom region of the image is zero, those pixels do not contribute to the inner products in the projections, and so the basis is effectively orthonormal.  Since the basis vectors are *not* zero in the atom region, when we take their sum weighted by their projections, we get the full reconstructed background image even in the atom region.  This is then used with the usual expression OD=log( abs( back./raw ) ) to get the optical depth.

I'm sure I'm not the first person to run into this problem so there are probably technical terms to describe what is done here, either from linear algebra or from image processing.  However I'm not an expert in either of those fields.  I have the feeling that this process can be desribed by a considering setting pixels to zero as acting on the vector with some metric, but I'm not really sure if that's a helpful or even accurate way to see this.  Feel free to email me if you know more and have a better way to interpret this, or maybe edit this file yourself and send a pull request.



# TODO:
- [ ] I think load\_image() assumes the image is 121-by-121 right now.  It shouldn't be hard to generalize it.
- [ ] Make process\_series\_\*() more general.  Maybe let it generate a file list from a filename\_pattern
- [ ] Add process\_series\_eig() to quick start section
- [ ] Prove eigenfaces are the eigenvectors with the largest eigenvalues

## To Maybe Do:
- [ ] Modify load\_image to call a user-supplied function to read the hard drive.  Expect the returned image to have the noise image subtracted
- [ ] Make object-oriented image processor class, so that the constructed basis, image size, etc. don't need to be passed around to each function.
- [ ] Ditch the svd algorithm completely
- [ ] Include a sample data set for testing.  Update docstring examples to use the sample data.
- [ ] Include a livescript with the examples from the doc strings
- [ ] Include notes from trying out different things in this Readme or some other file


# References

Xiaolin Li, Min Ke, Bo Yan, and Yuzhu Wang, "Reduction of interference fringes in absorption imaging of cold atom cloud using eigenface method," Chin. Opt. Lett. 5, 128-130 (2007)

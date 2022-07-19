import tensorflow as tf
# import tensorflow.compat.v1 as tf
import numpy as np

flags = tf.app.flags
FLAGS = flags.FLAGS
# -*- coding: utf-8 -*-

class DecagonOptimizer(object):
    def __init__(self, embeddings, latent_inters, latent_varies,
                 degrees, edge_types, edge_type2dim, placeholders,
                 margin=0.1, neg_sample_weights=1., batch_size=100):
        self.embeddings= embeddings
        self.latent_inters = latent_inters
        self.latent_varies = latent_varies
        self.edge_types = edge_types
        self.degrees = degrees
        self.edge_type2dim = edge_type2dim
        obj_type2n = {i: edge_type2dim[i,j][0][0] for i, j in edge_types}
        self.margin = margin
        self.neg_sample_weights = neg_sample_weights
        self.batch_size = batch_size

        inputs = placeholders['batch']
        batch_edge_type_idx = placeholders['batch_edge_type_idx']
        batch_row_edge_type = placeholders['batch_row_edge_type']
        batch_col_edge_type = placeholders['batch_col_edge_type']
        row_inputs = tf.squeeze(gather_cols(inputs, [0]))
        col_inputs = tf.squeeze(gather_cols(inputs, [1]))

        obj_type_n = [obj_type2n[i] for i in range(len(embeddings))]
        obj_type_lookup_start = tf.cumsum([0] + obj_type_n[:-1])
        obj_type_lookup_end = tf.cumsum(obj_type_n)
        labels = tf.reshape(tf.cast(row_inputs, dtype=tf.int64), [batch_size, 1])
        print('-----------------------------------------------')
        neg_samples_list = []
        for i, j in self.edge_types:
            for k in range(self.edge_types[i,j]):
                #if self.edge_types == (1,0):
                print('i', i)
                print('j', j )
                print('k,', k)
                print('k: ', k)
                neg_samples, _, _ = tf.nn.fixed_unigram_candidate_sampler(
                    true_classes=labels,
                    num_true=1,
                    num_sampled=batch_size,
                    unique=False,
                    range_max=len(degrees[0][i][k]),
                    distortion=0.75,
                    unigrams=degrees[0][i][k].tolist())
                print('neg_samples item: ', neg_samples)
                print('neg_samples device', neg_samples.device)
                print('neg_samples shape', neg_samples.get_shape)
                #print('neg samples eval: ', neg_samples.eval())
                # with tf.Session() as sess:     
                #     x_value = sess.run(neg_samples)
                #     print(neg_samples) #[1.  1.5 2. ]
                #     np.save("x.npy", neg_samples, allow_pickle=False)
                #one_string = tf.strings.format("{neg_samples}\n", (neg_samples))
                #print('test saving')
                #tf.io.write_file('test_saving_test.npy', one_string, name=None)
                #break
                # sess = tf.Session()
                # with sess.as_default():
                #     print(type(tf.constant([1,2,3]).eval()))
                #     print(type(neg_samples.eval()))
                #     print(neg_samples.eval())
                #sess = tf.InteractiveSession()
                # print(type(tf.constant([1,2,3]).eval()))
                print('...')
                #print('neg_samples to array', dir(neg_samples))

                neg_samples_list.append(neg_samples)
        print('neg_samples_list: ', neg_samples_list)
        print('len neg samples: ', len(neg_samples_list))
        print('ended ..........')
        
        self.neg_samples = tf.gather(neg_samples_list, self.batch_edge_type_idx)
        # aqui predice los que en teoria son positivos
        self.preds = self.batch_predict(self.row_inputs, self.col_inputs)
        self.outputs = tf.diag_part(self.preds)
        self.outputs = tf.reshape(self.outputs, [-1])
        # aqui los k en teoria son negativos
        
        self.neg_preds = self.batch_predict(self.neg_samples, self.col_inputs)
        self.neg_outputs = tf.diag_part(self.neg_preds)
        self.neg_outputs = tf.reshape(self.neg_outputs, [-1])

        self.predict()

        self._build()

    def batch_predict(self, row_inputs, col_inputs):
        concatenated = tf.concat(self.embeddings, 0)

        ind_start = tf.gather(self.obj_type_lookup_start, self.batch_row_edge_type)
        ind_end = tf.gather(self.obj_type_lookup_end, self.batch_row_edge_type)
        indices = tf.range(ind_start, ind_end)
        row_embeds = tf.gather(concatenated, indices)
        row_embeds = tf.gather(row_embeds, row_inputs)

        ind_start = tf.gather(self.obj_type_lookup_start, self.batch_col_edge_type)
        ind_end = tf.gather(self.obj_type_lookup_end, self.batch_col_edge_type)
        indices = tf.range(ind_start, ind_end)
        col_embeds = tf.gather(concatenated, indices)
        col_embeds = tf.gather(col_embeds, col_inputs)

        latent_inter = tf.gather(self.latent_inters, self.batch_edge_type_idx)
        latent_var = tf.gather(self.latent_varies, self.batch_edge_type_idx)

        product1 = tf.matmul(row_embeds, latent_var)
        product2 = tf.matmul(product1, latent_inter)
        product3 = tf.matmul(product2, latent_var)
        preds = tf.matmul(product3, tf.transpose(col_embeds))
        return preds

    def predict(self):
        concatenated = tf.concat(self.embeddings, 0)

        ind_start = tf.gather(self.obj_type_lookup_start, self.batch_row_edge_type)
        ind_end = tf.gather(self.obj_type_lookup_end, self.batch_row_edge_type)
        indices = tf.range(ind_start, ind_end)
        row_embeds = tf.gather(concatenated, indices)

        ind_start = tf.gather(self.obj_type_lookup_start, self.batch_col_edge_type)
        ind_end = tf.gather(self.obj_type_lookup_end, self.batch_col_edge_type)
        indices = tf.range(ind_start, ind_end)
        col_embeds = tf.gather(concatenated, indices)


        latent_inter = tf.gather(self.latent_inters, self.batch_edge_type_idx)
        latent_var = tf.gather(self.latent_varies, self.batch_edge_type_idx)

        product1 = tf.matmul(row_embeds, latent_var)
        product2 = tf.matmul(product1, latent_inter)
        product3 = tf.matmul(product2, latent_var)
        self.predictions = tf.matmul(product3, tf.transpose(col_embeds))

    def _build(self):
        # self.cost = self._hinge_loss(self.outputs, self.neg_outputs)
        self.cost = self._xent_loss(self.outputs, self.neg_outputs)
        self.optimizer = tf.train.AdamOptimizer(learning_rate=FLAGS.learning_rate)

        self.opt_op = self.optimizer.minimize(self.cost)
        self.grads_vars = self.optimizer.compute_gradients(self.cost)

    def _hinge_loss(self, aff, neg_aff):
        """Maximum-margin optimization using the hinge loss."""
        diff = tf.nn.relu(tf.subtract(neg_aff, tf.expand_dims(aff, 0) - self.margin), name='diff')
        loss = tf.reduce_sum(diff)
        return loss

    def _xent_loss(self, aff, neg_aff):
        """Cross-entropy optimization."""

        # l2_loss = tf.add_n(tf.get_collection("l2_reg"))
        # tf.add_to_collection('l2_reg', tf.contrib.layers.l2_regularizer(1.0)(self.vars['weights_%d' % k]))
        true_xent = tf.nn.sigmoid_cross_entropy_with_logits(labels=tf.ones_like(aff), logits=aff)
        negative_xent = tf.nn.sigmoid_cross_entropy_with_logits(labels=tf.zeros_like(neg_aff), logits=neg_aff)
        loss = tf.reduce_sum(true_xent) + self.neg_sample_weights * tf.reduce_sum(negative_xent)
        return loss


def gather_cols(params, indices, name=None):
    """Gather columns of a 2D tensor.

    Args:
        params: A 2D tensor.
        indices: A 1D tensor. Must be one of the following types: ``int32``, ``int64``.
        name: A name for the operation (optional).

    Returns:
        A 2D Tensor. Has the same type as ``params``.
    """
    with tf.op_scope([params, indices], name, "gather_cols") as scope:
        # Check input
        params = tf.convert_to_tensor(params, name="params")
        indices = tf.convert_to_tensor(indices, name="indices")
        try:
            params.get_shape().assert_has_rank(2)
        except ValueError:
            raise ValueError('\'params\' must be 2D.')
        try:
            indices.get_shape().assert_has_rank(1)
        except ValueError:
            raise ValueError('\'params\' must be 1D.')

        # Define op
        p_shape = tf.shape(params)
        p_flat = tf.reshape(params, [-1])
        i_flat = tf.reshape(tf.reshape(tf.range(0, p_shape[0]) * p_shape[1],
                                       [-1, 1]) + indices, [-1])
        return tf.reshape(
            tf.gather(p_flat, i_flat), [p_shape[0], -1])


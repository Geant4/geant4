"""
** model **
defines the VAE model class 
"""

# Setup
import keras
from tensorflow.keras.layers import Input, Dense, Lambda, Layer, Multiply, Add, concatenate
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.models import Model
from tensorflow.keras import backend as K
from tensorflow.keras import metrics

# VAE model class
class VAE:
    def __init__(self, **kwargs):
        self.original_dim = kwargs.get('original_dim')
        self.latent_dim = kwargs.get('latent_dim')
        self.batch_size = kwargs.get('batch_size')
        self.intermediate_dim1 = kwargs.get('intermediate_dim1')
        self.intermediate_dim2 = kwargs.get('intermediate_dim2')
        self.intermediate_dim3 = kwargs.get('intermediate_dim3')
        self.intermediate_dim4 = kwargs.get('intermediate_dim4')
        self.epsilon_std = kwargs.get('epsilon_std')
        self.mu = kwargs.get('mu')
        self.lr = kwargs.get('lr')
        self.epochs = kwargs.get('epochs')
        self.activ = kwargs.get('activ')
        self.outActiv = kwargs.get('outActiv')
        self.validation_split = kwargs.get('validation_split')
        self.wReco = kwargs.get('wReco')
        self.wkl = kwargs.get('wkl')
        self.optimizer = kwargs.get('optimizer')
        self.ki = kwargs.get('ki')
        self.bi = kwargs.get('bi')
        self.checkpoint_dir = kwargs.get('checkpoint_dir')
        self.earlyStop = kwargs.get('earlyStop')
        # KL divergence computation
        class KLDivergenceLayer(Layer):
            def __init__(self, *args, **kwargs):
                self.is_placeholder = True
                super(KLDivergenceLayer, self).__init__(*args, **kwargs)
            def call(self, inputs):
                mu, log_var = inputs
                kl_batch = -self.wkl * K.sum(1 + log_var - K.square(mu) - K.exp(log_var), axis=-1)
                self.add_loss(K.mean(kl_batch), inputs=inputs)
                return inputs
        # Build the encoder
        xIn = Input((input_dim,))
        eCond = Input(shape=(1,))
        angleCond = Input(shape=(1,))
        GeoCond = Input(shape=(2,))
        mergedInput = concatenate([xIn, eCond, angleCond, GeoCond],)
        h1 = Dense(self.intermediate_dim1, activation=self.activ,
                   kernel_initializer=self.ki, bias_initializer=self.bi)(mergedInput)
        h1 = BatchNormalization()(h1)
        h2 = Dense(self.intermediate_dim2, activation=self.activ,
                   kernel_initializer=self.ki, bias_initializer=self.bi)(h1)
        h2 = BatchNormalization()(h2)
        h3 = Dense(self.intermediate_dim3, activation=self.activ,
                   kernel_initializer=self.ki, bias_initializer=self.bi)(h2)
        h3 = BatchNormalization()(h3)
        h4 = Dense(self.intermediate_dim4, activation=self.activ,
                   kernel_initializer=self.ki, bias_initializer=self.bi)(h3)
        h = BatchNormalization()(h4)
        z_mu = Dense(self.latent_dim,)(h)
        z_log_var = Dense(self.latent_dim,)(h)
        # compute the KL divergence
        z_mu, z_log_var = KLDivergenceLayer()([z_mu, z_log_var])
        # Reparameterization trick
        z_sigma = Lambda(lambda t: K.exp(.5*t))(z_log_var)
        eps = Input(tensor=K.random_normal(shape=(K.shape(xIn)[0], self.latent_dim)))
        z_eps = Multiply()([z_sigma, eps])
        z = Add()([z_mu, z_eps])
        zCond = concatenate([z,eCond,angleCond,GeoCond],)
        # This defines the encoder which takes noise and input and outputs the latent variable z
        self.encoder = Model(inputs=[xIn,eCond,angleCond,GeoCond,eps], outputs=zCond)
        # Build the decoder / Generator
        decoL4 = Dense(self.intermediate_dim4, input_dim=(self.latent_dim+4),
                       activation=self.activ, kernel_initializer=self.ki, bias_initializer=self.bi)
        decoL4_BN = BatchNormalization()
        decoL3 = Dense(self.intermediate_dim3, input_dim=self.intermediate_dim4,
                       activation=self.activ, kernel_initializer=self.ki, bias_initializer=self.bi)
        decoL3_BN = BatchNormalization()
        decoL2 = Dense(self.intermediate_dim2, input_dim=self.intermediate_dim3,
                       activation=self.activ, kernel_initializer=self.ki, bias_initializer=self.bi)
        decoL2_BN = BatchNormalization()
        decoL1 = Dense(self.intermediate_dim1, input_dim=self.intermediate_dim2,
                       activation=self.activ, kernel_initializer=self.ki, bias_initializer=self.bi)
        decoL1_BN = BatchNormalization()
        x_reco = Dense(self.original_dim, activation=self.outActiv)
        zDecoInput = Input(shape=(latent_dim+4,))
        x_recoDeco = x_reco((((decoL1_BN(decoL1(decoL2_BN(decoL2(decoL3_BN(decoL3(decoL4_BN(decoL4(zDecoInput))))))))))))
        # This defines the decoder which takes an input of size latent dimension + condition size dimension and outputs the  reconstructed input version
        self.decoder = Model(inputs=[zDecoInput], outputs=[x_recoDeco]) 
        # This defines the reconstruction loss of the VAE model
        def reconstructionLoss(G4_Event, VAE_Event):
            return K.mean(self.wReco*K.sum(metrics.binary_crossentropy(G4_Event, VAE_Event)))
        # This defines the VAE model (encoder and decoder)
        self.vae = Model(inputs=[xIn,eCond,angleCond,GeoCond,eps], outputs=[self.decoder(self.encoder([xIn, eCond,angleCond,GeoCond,eps]))])
        self.vae.compile(optimizer=self.optimizer, loss=[reconstructionLoss] )
    # Training function 
    def train(self, trainSet, eCond, angleCond, GeoCond):
        # If the early stopping flag is on then stop the training when a monitored metric (validation) has stopped improving after (patience) number of epochs  
        if(self.earlyStop):
            from tensorflow.keras.callbacks import EarlyStopping
            cP = EarlyStopping(monitor='val_loss', min_delta=0.01, patience=5,verbose=1)
        # If the early stopping flag is off then run the training for the number of epochs and save the model every (period) epochs
        else:
            cP = keras.callbacks.ModelCheckpoint('%s/VAE-{epoch:02d}.h5'%self.checkpoint_dir, monitor='val_loss',
                                                verbose=0, save_best_only=False, save_weights_only=False, mode='auto',
                                                period=100)
        noise = np.random.normal(0,1, size = (trainSet.shape[0],latent_dim))
        history = self.vae.fit([trainSet, eCond, angleCond, GeoCond,noise], [trainSet], 
                               shuffle=True,
                               epochs=self.epochs,
                               verbose=1,
                               validation_split=self.validation_split,
                               batch_size=self.batch_size,
                               callbacks=[cP]
                               )
        return history
    # Encode function uses only the encoder to generate the latent representation of an input
    def encode(self, dataSet):
        return self.encoder.predict(dataSet, batch_size=self.batch_size)
    # Generate function uses only the decoder to generate new showers using the z_sample which is a vector of 10D Gaussians in addition to 
    def generate(self, z_sample):
        return self.decoder.predict([z_sample])
    # Encode function
    def predict(self, dataSet):
        return self.vae.predict(dataSet, batch_size=self.batch_size)
    # Encode function
    def evaluate(self, dataSet):
        return self.vae.evaluate(dataSet, batch_size=self.batch_size)

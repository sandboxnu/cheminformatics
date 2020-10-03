from flask import Flask
from flask_bootstrap import Bootstrap
import flask_excel as excel
from flask_session import Session

app = Flask(__name__, static_url_path='', static_folder='static')
app.config['PRODUCTION'] = False

bootstrap = Bootstrap()
bootstrap.init_app(app)

app.config['SESSION_TYPE'] = 'filesystem'
app.config['SECRET_KEY'] = 'some secret key' 
sess = Session()
sess.init_app(app)

excel.init_excel(app)

from app import routes
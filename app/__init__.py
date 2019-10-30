from flask import Flask
from flask_bootstrap import Bootstrap

app = Flask(__name__, static_url_path='', static_folder='static')

bootstrap = Bootstrap()
bootstrap.init_app(app)

from app import routes
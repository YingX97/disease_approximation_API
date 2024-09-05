from config import configuration as config

current_version = config['api_version']


def get_api_endpoint(resource_name, api_version='latest'):
    """Get the endpoint for an API resource"""
    if api_version == 'latest':
        api_version = current_version

    return f'/{api_version}/{resource_name}'